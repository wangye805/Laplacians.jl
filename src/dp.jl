#=

distributed approxChol Laplacian solver by Ye Wang 2020.

this solver is based on Prof Spielman's single-threaded approximate cholesky

=#
import Metis;
#apply metis graph separator to do nest dissection
#input:
#   adj: graph adjacency matrix
#   numHier: number of hierarchy's dissected. The whole graph will be cutted into 2^numHier parts
#output:
#   P: permutation vector
#   numPorts: last section in the permutation vector are boundary points
function nestDissection(adj::SparseArrays.SparseMatrixCSC{Tval,Tind}, numHier::Tind) where {Tval, Tind}
    #a: hierarchical adj matrix of sub matrices
    a = Array{Array{SparseArrays.SparseMatrixCSC{Tval, Tind}, 1}, 1}(undef, numHier);
    #s: set of nodes indices in global id
    s = Array{Array{Array{Tind}, 1}, 1}(undef, numHier);
    #b: set of boundary nodes
    b = Tind[];

    a[1] = [adj];
    s[1] = [collect(1:adj.n)];
    for ii=1:1:(numHier - 1)
        #at ii-th hier, we have 2^(ii-1) sub-matrix to dissection
        #prepare a[ii+1]
        a[ii+1] = Array{SparseArrays.SparseMatrixCSC{Tval, Tind}, 1}(undef, 2^ii);
        s[ii+1] = Array{Array{Tind}, 1}(undef, 2^ii);
        for jj=1:1:2^(ii-1)
            part = Metis.separator(a[ii][jj]);
            #find local indices of part1, part2 and boundry parts
            sLocal1 = findall(part.==1);
            sLocal2 = findall(part.==2);
            sLocal3 = findall(part.==3);
            a[ii+1][2*jj-1] = Laplacians.submatrix(a[ii][jj], sLocal1);
            s[ii+1][2*jj-1] = (s[ii][jj])[sLocal1];
            a[ii+1][2*jj] = Laplacians.submatrix(a[ii][jj], sLocal2);
            s[ii+1][2*jj] = (s[ii][jj])[sLocal2];
            b = [(s[ii][jj])[sLocal3]; b];
        end
    end
    #collect permutation matrix at the last hierachy
    P = Tind[];
    for ii = 1:1:2^(numHier-1)
        P = [P; s[numHier][ii]];
    end
    #check nodes are indeed part into 2^numHier parts
    adjP = Laplacians.submatrix(adj, P)
    @assert(maximum(Laplacians.components(adjP))>=2^(numHier-1));
    #now include ports
    P = [P; b];
    #check whether P is a perm
    @assert(isperm(P));
    return P, length(b);
end

#forwardSolve with multi partitions, used in partitioned version of solve!
function forwardSolve!(ldls::Array{LDL{Tind, Tval}, 1}, #result ldl obtained from approximate factorization
                       y::Vector, #the vector to solve
                       numInternalNodes::Array{Tind,1}, #num of internal nodes in each part
                       indexOffsets::Array{Tind, 1}, #index offset for each partition
                       portVecs::Array{Array{Tind,1}, 1}, #port mapping from each part to aggregated ports
                       numInternalNode::Tind,
                       numPort::Tind)where {Tind, Tval}
    #get num of partitions
    numPart = size(ldls, 1);
    #check partition size
    @assert(numPart == size(portVecs, 1));
    @assert(numPart == size(indexOffsets, 1));
    
    for ii in 1:numPart
        #get array slicing for each partition
        yiInternal = y[indexOffsets[ii]+1:indexOffsets[ii]+numInternalNodes[ii]];
        #check port dimension
        @assert(ldls[ii].m==size(portVecs[ii],1));
        yiPort = y[portVecs[ii].+numInternalNode];
        yi = [yiInternal;yiPort];
        #call forward solve in each partition
        forwardSolve!(ldls[ii], yi); 
        #assign yi back to the aggregated vector y
        y[indexOffsets[ii]+1:indexOffsets[ii]+numInternalNodes[ii]] = yi[1:numInternalNodes[ii]];
        y[portVecs[ii].+numInternalNode] = yi[numInternalNodes[ii]+1:end];
    end
    return y;
end

#backwardSolve with multi partitions, used in partitioned version of solve
function backwardSolve!(ldls::Array{LDL{Tind, Tval}, 1}, #result ldl obtained from approximate factorization
                       y::Vector, #the vector to solve
                       numInternalNodes::Array{Tind,1}, #num of internal nodes in each part
                       indexOffsets::Array{Tind, 1}, #index offset for each partition
                       portVecs::Array{Array{Tind,1}, 1}, #port mapping from each part to aggregated ports
                       numInternalNode::Tind,
                       numPort::Tind)where {Tind, Tval}
    #get num of partitions
    numPart = size(ldls, 1);
    #check partition size
    @assert(numPart == size(portVecs, 1));
    @assert(numPart == size(indexOffsets, 1));
    
    for ii in 1:numPart
        #get array slicing for each partition
        yiInternal = y[indexOffsets[ii]+1:indexOffsets[ii]+numInternalNodes[ii]];
        #check port dimension
        @assert(ldls[ii].m==size(portVecs[ii],1));
        yiPort = y[portVecs[ii].+numInternalNode];
        yi = [yiInternal;yiPort];
        #call backward solve in each partition
        backwardSolve!(ldls[ii], yi); 
        #assign yi back to the aggregated vector y
        y[indexOffsets[ii]+1:indexOffsets[ii]+numInternalNodes[ii]] = yi[1:numInternalNodes[ii]];
        y[portVecs[ii].+numInternalNode] = yi[numInternalNodes[ii]+1:end];
    end
    return y;
end

"""
Solve (LDL)^-1 with approximate factorization results obtained
multiple partition version
input: 
    ldls: array of LDL, approximate factorization results
    ls: graph laplacian of the aggregated schur complement
    portVecs: mapping from ports in each partition to aggregated ports
    numPort: total number of ports in the entire system
input&output:
    y: vector for solving
"""
function solve!(ldls::Array{LDL{Tind, Tval}, 1}, 
                y::Vector, 
                ls::SparseMatrixCSC{Tval, Tind}, 
                numInternalNodes::Array{Tind,1}, #num of internal nodes in each part
                indexOffsets::Array{Tind,1}, 
                portVecs::Array{Array{Tind,1}, 1}, 
                numInternalNode::Tind, 
                numPort::Tind)where {Tind, Tval}
    forwardSolve!(ldls, y, numInternalNodes, indexOffsets, portVecs, numInternalNode, numPort);
    #check number of internal nodes
    @assert((numInternalNode + numPort) == size(y, 1));
    #schur^{-1}
    #check schur laplacian size
    @assert(size(ls,1)==numPort);
    S = ls[1:end-1, 1:end-1];
    if numPort>1
        x2 = [S\y[numInternalNode+1:numInternalNode + numPort-1]; 0];
    else
        x2 = [0];
    end
    y[numInternalNode+1:numInternalNode + numPort] = x2;

    backwardSolve!(ldls, y, numInternalNodes, indexOffsets, portVecs, numInternalNode, numPort);
    return y;
end

"""
The function to do laplacian matrix vector multiplication y = LA*v with partitioned format
"""
function LaplacianVectorMult!(las::Array{SparseMatrixCSC{Tval, Tind}, 1}, #graph laplacians for adj graph of each partition, here we include the gnd node at last 
                               y::Vector, #input vector for multiplication
                               numInternalNodes::Array{Tind,1}, #num of internal nodes in each part
                               indexOffsets::Array{Tind,1}, #starting index of each partition
                               portVecs::Array{Array{Tind,1}, 1}, #port mapping for each partitions
                               numInternalNode::Tind, #total number of internal nodes
                               numPort::Tind)where {Tind, Tval}
    #get num of partitions
    numPart = size(las, 1);
    #check partition size
    @assert(numPart == size(portVecs, 1));
    #check vector y size
    @assert(size(y,1)==(numInternalNode + numPort));
    #create empty result so that later we can do mac operations
    x = copy(y);
    #y = zeros(size(x));
    fill!(y, 0.0);

    for ii in 1:numPart
        xiInternal = x[indexOffsets[ii]+1:indexOffsets[ii]+numInternalNodes[ii]];
        xiPort = x[portVecs[ii].+numInternalNode];
        xi = [xiInternal;xiPort];
        yi = las[ii]*xi;
        y[indexOffsets[ii]+1:indexOffsets[ii]+numInternalNodes[ii]] = y[indexOffsets[ii]+1:indexOffsets[ii]+numInternalNodes[ii]] + yi[1:numInternalNodes[ii]];
        y[portVecs[ii].+numInternalNode] = y[portVecs[ii].+numInternalNode] + yi[numInternalNodes[ii]+1:end];
    end
    return y;
end

"""
The function to get the full adj graph from partitioned version
"""
function fullAdjGraph(adjGraphs::Array{SparseMatrixCSC{Tval, Tind}, 1}, #adj graph for each partition
                  portVecs::Array{Array{Tind, 1}, 1}, #port mapping vectors (from each partition to aggregated ports)
                  numPort::Tind; #total number of ports
                  verbose=false) where {Tind, Tval}

    #check partition size
    numPart = size(adjGraphs, 1);
    @assert(numPart == size(portVecs, 1));
    #first calculate problem size: total # of nodes, # of ports, # of internal nodes
    numNodes = Array{Tind}(undef, numPart);
    numPorts = Array{Tind}(undef, numPart);
    for ii in 1:numPart
        numNodes[ii] = adjGraphs[ii].m;#adjGraphs[ii].m should be equal to adjGraphs[ii].n
        numPorts[ii] = size(portVecs[ii], 1);       
    end
    numInternalNodes = numNodes .- numPorts;
    indexOffsets = accumulate(+, [0; numInternalNodes[1:end-1]]);
    numInternalNode = sum(numInternalNodes);

    I = Tind[];
    J = Tind[];
    V = Tind[];
    #initialize an empty result adjGraph
    for ii in 1:numPart
        (Iii, Jii, Vii) = findnz(adjGraphs[ii]);
        @inbounds for jj in 1:size(Vii, 1)
            #process I vector
            if Iii[jj]>numInternalNodes[ii] #port node
                Iii[jj] = numInternalNode + (portVecs[ii])[Iii[jj] - numInternalNodes[ii]];
            else #internal node
                Iii[jj] = Iii[jj] + indexOffsets[ii];
            end
            #same for J vector
            if Jii[jj]>numInternalNodes[ii] #port node
                Jii[jj] = numInternalNode + (portVecs[ii])[Jii[jj] - numInternalNodes[ii]];
            else #internal node
                Jii[jj] = Jii[jj] + indexOffsets[ii];
            end
        end
        I = [I; Iii];
        J = [J; Jii];
        V = [V; Vii];
    end
    adjGraph = sparse(I, J, V, numInternalNode+numPort, numInternalNode+numPort)
    return adjGraph;
end

"""
The function to aggregate partioned schur complements to a whole one
"""
function schurComplement(schurCs::Array{SparseMatrixCSC{Tval, Tind}, 1}, #schur complement obtained through approximate factorization in each partition
                         portVecs::Array{Array{Tind, 1}, 1}, #port mapping vectors (from each partition to aggregated ports
                         numPort::Tind) where {Tind, Tval}
    #check partition size
    numPart = size(schurCs, 1);
    @assert(numPart == size(portVecs, 1));
    I = Tind[];
    J = Tind[];
    V = Tind[];
    #initialize an empty result schur complement matrix
    schurC = sparse(Tind[], Tind[], Tval[], numPort, numPort)
    for ii in 1:numPart
        (Iii, Jii, Vii) = findnz(schurCs[ii]);
        I = Array{Tind}(undef, size(Iii, 1));
        J = Array{Tind}(undef, size(Jii, 1));
        #now map correct index from Iii, Jii to I, J
        I = (portVecs[ii])[Iii];
        J = (portVecs[ii])[Jii];
        subSchurC = sparse(I, J, Vii, numPort, numPort);
        schurC = schurC + subSchurC;
    end
    return schurC;
end


"""
condition number for partitioned approximate factorization
"""
function condNumber(adjGraphs::Array{SparseMatrixCSC{Tval, Tind}, 1}, #adj graph for each partition
                    ldls::Array{LDL{Tind, Tval}, 1}, #LDL results obtained through approximate factorization in each partition 
                    schurCs::Array{SparseMatrixCSC{Tval, Tind}, 1}, #schur complement obtained through approximate factorization in each partition
                    portVecs::Array{Array{Tind, 1}, 1}, #port mapping vectors (from each partition to aggregated ports)
                    numPort::Tind; #total number of ports
                    verbose=false) where {Tind, Tval}
    #check partition size
    numPart = size(adjGraphs, 1);
    @assert(numPart == size(ldls, 1));
    @assert(numPart == size(schurCs, 1));
    @assert(numPart == size(portVecs, 1));

    #first calculate problem size: total # of nodes, # of ports, # of internal nodes
    numNodes = Array{Tind}(undef, numPart);
    numPorts = Array{Tind}(undef, numPart);
    #create array of laplacian matrices
    las = Array{SparseMatrixCSC{Tval, Tind}}(undef, numPart);
    for ii in 1:numPart
        numNodes[ii] = adjGraphs[ii].m;#adjGraphs[ii].m should be equal to adjGraphs[ii].n
        numPorts[ii] = size(portVecs[ii], 1);       
        #create graph laplacian for each partition
        las[ii] = lap(adjGraphs[ii]);
    end
    numInternalNodes = numNodes .- numPorts;
    indexOffsets = accumulate(+, [0; numInternalNodes[1:end-1]]);
    numInternalNode = sum(numInternalNodes);
    #also create the aggreated schur complement and its corresponding laplacian
    schurC = schurComplement(schurCs, portVecs, numPort);
    ls = lap(schurC);

    #create the square linear operation
    g = function(b)
        y = copy(b);
        #first do the matrix vector multiplication but with many partitions
        y = LaplacianVectorMult!(las, y, numInternalNodes, indexOffsets, portVecs, numInternalNode, numPort);
        #then do the forward solve with given factorization results
        solve!(ldls, y, ls, numInternalNodes, indexOffsets, portVecs, numInternalNode, numPort);
        return y;
    end

    #total number of node is sum of internal nodes + port nodes
    numNode  = numInternalNode + numPort;
    gOp = SqLinOp(false,1.0,numNode,g)
    upper = eigs(gOp;nev=1,which=:LM,tol=1e-2)[1][1]
    
    g2(b) = abs(upper)*b - g(b)
    g2Op = SqLinOp(false,1.0,numNode,g2)
    lower = abs(upper) - eigs(g2Op;nev=2,which=:LM,tol=1e-2)[1][2]
    
    if verbose
        println("lower: ", lower, ", upper: ", upper);
    end
    
    return upper/lower
end

