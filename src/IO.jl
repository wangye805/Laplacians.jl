"""
    adj = read_graph(fn)

Read a graph from a file in IJ or IJV format.
That is, each line of the file should represent one edge.
Each edge should be specified by the indices of its vertices,
separated by a whitespace or a comma.  If the graph is weighted,
the weight should follow the second index.  For example, the unweighted
complete graph on 3 vertices would appear as the file

```
1 2
1 3
2 3
```

A weighted path on 3 vertices with edge weights `1.5` and `2.5` would be

```
1, 2, 1.5
2, 3, 2.5
```

The function tries to detect the delimiter type (comma or whitespace)
from the first line of the file.  The format must be consistent.
Vertex indices start at `1`.

Input matrix is only upper or lower triangular part of the graph adj matrix
"""
function read_graph(fn::AbstractString)

    # first, detect the delimiter type
    fh = open(fn)
    ln = readline(fh)
    close(fh)

    if occursin(",", ln)
        data = readdlm(fn,',')
    else
        data = readdlm(fn)
    end



    n = maximum(data[:,1:2])

    edlist = convert(Array{Int64,2}, data[:,1:2])

    if size(data,2) == 3
        wts = convert(Array{Float64,1},data[:,3])
    else
        wts = ones(size(data,1))
    end

    a = sparse(edlist[:,1],edlist[:,2],wts,n,n)
    a = a + a'

    return a
end

#the function to read matrix market format in which:
#the first line specify the matrix (graph) size (m, n, nnz)
#then each line specify a matrix entry for (i>=j).
#further we assume diagonal entries proceed off-diagonal entries
#output of this function is a MNA matrix A (removing gnd node from the corresponding Laplacian)
function read_matrix_market(fn::AbstractString)
	data = readdlm(fn)
	#the first line specify the m, n, num of lines
	#matrix size is of mxn
	m = convert(Int64, data[1,1]);
	n = convert(Int64, data[1,2]);
	if m!=n
		println("matrix size m neq n");
	end
	#num of lines contains diagonal (i==j) + lower triangular part (i>j)
	num_lines = convert(Int64, data[1,3]);
	if(size(data,1)!=(num_lines+1))
		println("file line size inconsistent");
	end
	#extract the diagonal part
	diag_list = convert(Array{Int64, 2}, data[2:1+m, 1:2]);
	diag_wt = convert(Array{Float64, 1}, data[2:1+m, 3]);
	#extract the lower triangular part (i>j)
	lt_list = convert(Array{Int64, 2}, data[2+m:end, 1:2]);
	lt_wt = convert(Array{Float64, 1}, data[2+m:end, 3]);
	a = sparse(lt_list[:,1], lt_list[:,2], lt_wt, n, n);
	a = a + a';
	a = a + sparse(diag_list[:,1], diag_list[:,2], diag_wt, n, n);
end

#the function to read the rhs format from Lengfei
#basically we ignore the first column and sum up the entire second column to get gnd current value
function read_rhs(fn::AbstractString)
	data = readdlm(fn);
	#convert the second column into float 64
	b = convert(Array{Float64}, data[:,2]);
	return vcat(b, -sum(b));
end

#the function to translate matrix market to graph edge list
#each edge only appear once (lower triangular)
function matrix_market_to_graph(fn::AbstractString)
	#first we read in the file in matrix market format
	A = read_matrix_market(fn);
	if size(A,1) != size(A,2)
		println("can only deal with square matrices");
	end
	n = size(A, 1) + 1;
	#sum up rows of A to get the connection to the gnd node (reference node)
	hidden_wt = sum(A, dims=1);
	#set the threshold to 10^-9 for meaningful conductance
	hidden_wt[hidden_wt.<10^-9].=0;
	hidden_wt = sparse(hidden_wt);
	(hi, hj, hv) = findnz(hidden_wt);
	(ai, aj, av) = findnz(tril(A, -1));
	ai = vcat(ai, n*ones(Int64, size(hi,1)));
	aj = vcat(aj, hj);
	av = vcat(-av, hv);
	return ai, aj, av, n
end

"""Writes the upper portion of a matrix in ijv format, one row for each edge,
separated by commas.  Only writes the upper triangular portion.
The result can be read from Matlab like this:

```
>> dl = dlmread('graph.txt');
>> a = sparse(dl(:,1),dl(:,2),dl(:,3));
>> n = max(size(a))
>> a(n,n) = 0;
>> a = a + a';
```
"""
function writeIJV(filename::AbstractString, mat)

  (ai,aj,av) = findnz(triu(mat))
  fh = open(filename,"w")
  for i in 1:length(ai)
    write(fh, "$(ai[i]),$(aj[i]),$(av[i])\n")
  end
  close(fh)

end #write IJV

import SparseArrays;
#the function to write a Laplacian matrix as the matrix market format
function write_matrix_market(filename::AbstractString, mat::SparseMatrixCSC)
	n = size(mat,1);
	if n != size(mat,2)
		println("The given matrix dimensions are not consistent");
	end
	A = mat[1:end-1, 1:end-1];
	fh = open(filename, "w");
	#write the first line
	write(fh, "$(n-1) $(n-1) $(SparseArrays.nnz(triu(A)))\n");
	#write all diag elements
	diag_wt = diag(A, 0);#usually this is a dense vector (for connected graphs)
	diag_list = 1:n-1;
	for i in 1:n-1
		write(fh, "$(diag_list[i]) $(diag_list[i]) $(diag_wt[i])\n");
	end
	#now write the off diagonal elements
	(ai, aj, av) = findnz(tril(A, -1));
	for i in 1:length(ai)
		write(fh, "$(ai[i]) $(aj[i]) $(av[i])\n");
	end
	close(fh);
end
