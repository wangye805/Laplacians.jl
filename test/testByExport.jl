# Test every function exported in Laplacians.jl by running it at least once
# Create this by copying over the export commands from Laplacians.jl

n = 101
a = wted_chimera(n,1)

# export symPermuteCSC

rp = randperm(n)
ap = symPermuteCSC(a,rp);
rpi = zeros(Int,n);
rpi[rp] = 1:n

@test sum(abs.(a-symPermuteCSC(ap,rpi))) == 0

# export symTransposeCSC

a2 = triu(a) + (tril(a) .> 0)
@test sum(abs.(symTransposeCSC(a2) - a2')) == 0

# export submatrixCSC

s = collect(1:10)
submatrixCSC(a,s)

  # export deg
  # export nbri
  # export weighti
  # export nbrs
  # export wdeg

x = 0
y = 0
z = 0
w = 0
for i in 1:n
    for j in 1:deg(a,i)
        global x += weighti(a,i,j)
        global y += a[i,nbri(a,i,j)]
    end
    for j in nbrs(a,i)
        global z += a[i,j]
    end
    global w += wdeg(a,i)
end
@test isapprox(x,y)
@test isapprox(x,z)
@test isapprox(x,w)
@test isapprox(x,sum(a))

# export setValue

setValue(a,1,1,0.0)

# export backIndices

b = backIndices(a)

# export flipIndex

b = flipIndex(a)

# export findEntries

u,v,w = findEntries(a)

# export compConductance

compConductance(a,collect(1:10))

# export getVolume

getVolume(a,collect(1:10))

# export getObound

getObound(a,collect(1:10))

  # export readIJ
  # export ringGraph
  # export generalizedRing
  # export randMatching
  # export randRegular


# export grownGraph

a2 = grownGraph(100,3)

# export grownGraphD

a2 = grownGraphD(100,3)

# export prefAttach

a2 = prefAttach(100,3,0.5)
a2 = prefAttach(5,4,0.5)


# export hyperCube

a2 = hyperCube(3)

# export completeBinaryTree

a2 = completeBinaryTree(7)

# export completeGraph

a2 = completeGraph(7)

# export pathGraph

a2 = pathGraph(7)

# export grid2

a2 = grid2(3)

# export grid2coords

a2 = grid2coords(3)

# export grid3

a3 = grid3(3)

  # export randGenRing
  # export randperm

# export ErdosRenyi

a2 = ErdosRenyi(100,300)

# export ErdosRenyiCluster

a2 = ErdosRenyiCluster(100,4)

# export ErdosRenyiClusterFix

a2 = ErdosRenyiClusterFix(100,4)

# export pureRandomGraph

# export chimera
# export randWeight
# export wted_chimera, semiWtedChimera

for i in 1:5
    semiWtedChimera(10000,i)
end


# export readIJ, readIJV, writeIJV
println("Testing IO")
n = 101
a = wted_chimera(n,1)
writeIJV("tmp.txt",a)
a2 = readIJV("tmp.txt")
@test sum(abs.(a-a2)) == 0

a2 = read_graph("tmp.txt")
@test sum(abs.(a-a2)) == 0

fh = open("tmp.txt","w")
write(fh,"1, 3, 4 \n 2, 3, 2.5 \n")
close(fh)

a1 = read_graph("tmp.txt")

fh = open("tmp.txt","w")
write(fh,"1 3 4 \n 2 3 2.5 \n")
close(fh)

a2 = read_graph("tmp.txt")

@test a1 == a2

rm("tmp.txt")

#testing APIs related to matrix market format
#testcase in page 2 of icm10post.pdf by Daniel Spielman
Lio = [2.0 -1 0 0 -1; -1 3 -1 -1 0; 0 -1 2 -1 0; 0 -1 -1 4 -2; -1 0 0 -2 3];
LioSparse = SparseArrays.sparse(Lio);
write_matrix_market("tmp.txt", LioSparse);
#tmp.txt should be
#4 4 8
#1 1 2
#2 2 3
#3 3 2
#4 4 4
#2 1 -1
#3 2 -1
#4 2 -1
#4 3 -1

#now read back from tmp.txt
A = read_matrix_market("tmp.txt");
#the golden of A should be [2 -1 0 0; -1 3 -1 -1; 0 -1 2 -1; 0 -1 -1 4];
@test norm(A-Lio[1:4, 1:4],2)<=1e-10
#test writing graph to (i, j, v) format
(ai, aj, av, n) = matrix_market_to_graph("tmp.txt")
#golden of (ai, aj, av) should be
#1 2 1
#1 5 1
#2 3 1
#2 4 1
#3 4 1
#4 5 2

#create the adj matrix based on i,j,v
#n = 5;
Ar = sparse(ai, aj, av, n, n);
Ar = Ar + Ar';
@test norm(lap(Ar)-Lio, 2)<=1e-10
rm("tmp.txt")

Laplacians.write_tesla("tmp.txt", LioSparse);
#the golden should be
#1 1 2.0
#2 2 3.0
#3 3 2.0
#4 4 4.0
#1 2 -1.0
#2 3 -1.0
#2 4 -1.0
#3 4 -1.0

#now read back from tmp
A = Laplacians.read_tesla("tmp.txt", 4);
#the golden of A should be [2 -1 0 0; -1 3 -1 -1; 0 -1 2 -1; 0 -1 -1 4];
@test LinearAlgebra.norm(A-Lio[1:4, 1:4],2)<=1e-10
#test writing graph to (i, j, v) format
(ai, aj, av, n) = Laplacians.tesla_to_graph("tmp.txt", 4);
#golden of (ai, aj, av) should be
#2 1 1
#3 2 1
#4 2 1
#4 3 1
#5 1 1
#5 4 2

#create the adj matrix based on i,j,v
#n = 5;
Ar = SparseArrays.sparse(ai, aj, av, n, n);
Ar = Ar + Ar';
@test LinearAlgebra.norm(Laplacians.lap(Ar)-Lio, 2)<=1e-10
rm("tmp.txt")

#test permutation of submatrix
P = [4, 2, 1, 3, 5];
PA = Laplacians.submatrix(Ar, P);
IPA = [3, 2, 3, 5, 2, 1, 2, 4, 1, 4, 5, 1];
JPA = [2, 3, 5, 3, 1, 2, 4, 2, 4, 1, 1, 5];
VPA = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2.0];
@test PA==SparseArrays.sparse(IPA, JPA, VPA, 5, 5);

#test submatrix and permutation of submatrix
P = [4, 2, 1];
PA = Laplacians.submatrix(Ar, P);
IPA = [3, 2, 1, 2];
JPA = [2, 3, 2, 1];
VPA = [1, 1, 1, 1.0];
@test PA==SparseArrays.sparse(IPA, JPA, VPA, 3, 3);


# export unweight, unweight!

a2 = unweight(a2)
unweight!(a2)

  # export mapweight
  # export uniformWeight, uniformWeight!

a2 = uniformWeight(a2)
uniformWeight!(a2)

  # export edgeVertexMat

b = edgeVertexMat(a2)
@test sum(abs.(b'*b - lap(unweight(a2)))) == 0

a = wted_chimera(102,2)
b = wtedEdgeVertexMat(a)
@test sum(abs.(b'*b - lap(a))) < 1e-8

  # export power, thicken_once, thicken

  a = power(grid2(10),4)
  a = thicken_once(grid2(10))
  a = thicken(grid2(10),4)

  # export productGraph
  # export generalizedNecklace
  # export subsampleEdges

subsampleEdges(a,0.5)

# export twoLift

twoLift(a)
twoLift(a,3)

  # export joinGraphs, disjoin

  # export plotGraph

  # export shortIntGraph, floatGraph

  # export lap
  # export adj
  # export spectral_coords
  # export spectral_drawing

  # export toUnitVector

# export diagmat

diagmat(a)

  # export components
  # export biggestComp
  # export vecToComps
  # export isConnected

# export shortestPaths, shortestPathTree, pathFromParents

a = wted_chimera(102,1)
shortestPaths(a,1)
shortestPathTree(a,1)

Laplacians.intHeapSort(randn(10))

nh = Laplacians.intHeap(10)
for i in 1:10
    Laplacians.intHeapAdd!(nh, i, rand())
end
Laplacians.intHeapSort(nh)

  # export kruskal, prim


# export RootedTree
# export matToTree

# export matToTreeDepth

a = wted_chimera(101,1)
t = akpw(a)
tr = matToTree(t)
tr, d1 = matToTreeDepth(t);
d2 = Laplacians.treeDepthDFS(t)

  # export tarjanStretch
  # export compDepth
  # export comp_stretches
  # export dfsOrder

Laplacians.bfsOrder(t,1)

t0 = complete_binary_tree(6); t0[1,2] = 0; t0[2,1] = 0; dropzeros!(t0);
@test_throws ErrorException Laplacians.bfsOrder(t0, 1)
@test_throws ErrorException Laplacians.matToTree(t0)
@test_throws ErrorException Laplacians.matToTreeDepth(t0)

  # export cg, cgSolver
  # export pcg, pcgSolver, pcgLapSolver

  # export maxflow

for i in 1:10
  a = wted_chimera(100,i)
  a[90:99,100] .= 2;
  a[100,90:99] .= 2;
  a[1,2:10] .= 2;
  a[2:10,1] .= 2;
  f,c = maxflow(a,1,100)

  @test sum(abs.(f+f')) < 1e-8
  y = f*ones(100)
  y[1] = 0
  y[100] = 0
  @test sum(abs.(y)) < 1e-8

  x = zeros(100)
  x[c] .= 1.0
  t = findall(iszero,x)
  @test abs( sum(a[c,t]) - sum(f[1,:])  ) < 1e-8

  @test maximum(f - a) < 1e-8

end

  # export akpw, akpwU

akpw(wted_chimera(10000,1),ver=2)

  # export prn
  # export apr
  # export localImprove
  # export refineCut
  # export dumb

a = chimera(100, 3);
s = prn(a, [1,2,3], 0.2, 5);
conds = compConductance(a, s)
#println(conds, " ", length(s))
minEpsSigma = getVolume(a, s) / getVolume(a, setdiff(collect(1:max(a.n, a.m)), s));
cut, flow = localImprove(a, s, epsSigma = minEpsSigma);
condcut = compConductance(a, cut)
heur = refineCut(a, cut)
dumbRefineCut(a,collect(1:10))

#println(condcut, " ", length(cut))

  # export randishKruskal, randishPrim

  # export FastSampler, sample, sampleMany


r = rand(10)
s = FastSampler(r)
sample(s)
blockSample(r)

# export fiedler

fiedler(chimera(100))

  # export SolverTest, speedTestLapSolvers

solvers = [SolverTest(approxchol_lap,"ac") SolverTest(augTreeLap,"aug")]

dic = Dict()
n = 1000
a = chimera(n)
b = randn(n)
b = b .- mean(b)
x = speedTestLapSolvers(solvers, dic, a, b, tol=1e-2, verbose=true)

f = Laplacians.augTreeFactor(a, akpw(a));

  # include("conditionNumber.jl")
  # export support, approxQual, conditionNumber

  # include("sparsify.jl")
  # export sparsify

a = wted_chimera(1000,1)
dave = nnz(a)/size(a,1)
a = thicken(a,round(Int,200/dave))
as = sparsify(a,ep=1,JLfac=4);
@test approxQual(a,as,verbose=true) < 2
@test conditionNumber(a,as,tol=1e-4) < 10



  # export johnlind

a = chimera(n,1)
johnlind(a)

  # export toposort, dirEdgeVertexMat

a = wted_chimera(301,1)
S = [1;2;3]
vals = [1.0;2.0;3.0]
x = harmonic_interp(a, S, vals, tol=1e-10)
@test x[S] == vals
b = lap(a)*x
b[S] .= 0
@test sum(abs.(b)) < 1e-6

#test the early stop mode using Icm10Post
II = [1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5];
JJ = [2, 5, 1, 3, 4, 2, 4, 2, 3, 5, 1, 4];
VV = [1.0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2];

a = SparseArrays.sparse(II, JJ, VV);
llmat = Laplacians.LLmatp(a);
#numAddPorts is 1
ldl, schurC, P = Laplacians.approxChol(llmat, 2, 1);
ldlGolden = Laplacians.LDL{Int64, Float64}(2, 3, 4, [1, 2], 1, [1, 1, 1], Int64[], Float64[], 5, [1, 3, 5], [1, 2, 1, 3], [1.0, 1, 1, 1], [2.0, 2]);
@test Laplacians.isapprox(ldl, ldlGolden);
schurCGolden = SparseArrays.sparse([1, 2, 1, 3, 2, 3], [2, 1, 3, 1, 3, 2], [1.5, 1.5, 0.5, 0.5, 2, 2]);
@test isapprox(schurC, schurCGolden);
@test isapprox(P, [3, 1, 2, 4, 5]);
#test the condition number, should be 1
aP = a[P,P];
cn = Laplacians.condNumber(aP, ldl, schurC,verbose=true);
@test isapprox(abs(cn), 1);

#test the corner case of early stop
#numAddPorts = 0
llmat = Laplacians.LLmatp(a);
ldl, schurC, P = Laplacians.approxChol(llmat, 2, 0);
ldlGolden = Laplacians.LDL{Int64, Float64}(3, 2, 4, [1, 2, 3], 3, [1, 2, 3, 3], [3, 3], [1.0, 1.0], 5, [1, 2, 3, 5], [1, 2, 1, 2], [1.0, 1.0, 1.5, 0.5], [2.0, 2.0, 2.0]);
schurCGolden = SparseArrays.sparse([1, 2], [2, 1], [2.375, 2.375]);
@test Laplacians.isapprox(ldl, ldlGolden);
@test isapprox(schurC, schurCGolden);
@test isapprox(P, [3, 1, 2, 4, 5]);
aP = a[P,P];
cn = Laplacians.condNumber(aP, ldl, schurC,verbose=true);
@test isapprox(abs(cn), 1);


#no ports get reduced
llmat = Laplacians.LLmatp(a);
ldl, schurC, P = Laplacians.approxChol(llmat, 2, 3);
ldlGolden = Laplacians.LDL{Int64, Float64}(0, 5, 4, Int64[], 1, [1], Int64[], Float64[], 1, [1], Int64[], Float64[], Float64[]);
schurCGolden = SparseArrays.sparse([3, 4, 3, 5, 1, 2, 4, 1, 3, 5, 2, 4], [1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5], [1.0, 1, 1, 1, 1, 1,1,1,1,2,1,2]);
@test Laplacians.isapprox(ldl, ldlGolden);
@test isapprox(schurC, schurCGolden);
@test isapprox(P, [3, 1, 2, 4, 5]);
aP = a[P,P];
cn = Laplacians.condNumber(aP, ldl, schurC,verbose=true);
@test isapprox(abs(cn), 1);

#test the debug and fixed random flag
a = Laplacians.grid2(3);
llmat = Laplacians.LLmatp(a); 
#reference run (debug = true, fixed random = false)
ldl, schurC, debugInfo = Laplacians.approxChol(llmat, 2; debug = true)
@test abs(Laplacians.condNumber(a, ldl, schurC,verbose=true)) <=2
order = ldl.col;
reduceDegs, PRNGs = debugInfo;
#when debug = false, fixed random = true
llmat = Laplacians.LLmatp(a);
ldlFixed, schurCFixed, debugInfo = Laplacians.approxChol(llmat, 2; fixedRandom = true, order = order, PRNGs = PRNGs);
@test Laplacians.isapprox(ldl, ldlFixed);
@test isapprox(schurC, schurCFixed);
#when debug = true, fixed random = true
llmat = Laplacians.LLmatp(a);
ldlFixed, schurCFixed, debugInfo = Laplacians.approxChol(llmat, 2; debug = true, fixedRandom = true, order = order, PRNGs = PRNGs);
@test Laplacians.isapprox(ldl, ldlFixed);
@test isapprox(schurC, schurCFixed);

#test the symmetric support methods (more robust but not as accurate as condNumber)
lower = Laplacians.support(a, ldl, schurC);
upper = Laplacians.support(ldl, schurC, a);
@test lower>=1
@test upper>=1
cn = abs(Laplacians.condNumber(a, ldl, schurC,verbose=true));
@test abs(upper*lower-cn)/cn<0.05

#test restrict approx chol results in a subset, lets use icm10post
I0 = [1, 2, 1, 5, 2, 3, 2, 4, 3, 4, 4, 5];
J0 = [2, 1, 5, 1, 3, 2, 4, 2, 4, 3, 5, 4];
V0 = [1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2];
#now interleave two identical copies of icm10post
II = vcat(I0.*2, (I0.*2).-1);
JJ = vcat(J0.*2, (J0.*2).-1);
VV = vcat(V0, V0);
adjGraph = SparseArrays.sparse(II, JJ, VV);
#now test the approx chol with disconnected forests
llmat = Laplacians.LLmatp(adjGraph);
ldl, schurC = Laplacians.approxChol(llmat, 2);
subset = [1, 3, 5, 7, 9];
ldlSub, schurCSub = Laplacians.subsetApproxChol(ldl, schurC, subset);
#should have the same results with single copy
subAdjGraph = adjGraph[subset, subset];
llmat = Laplacians.LLmatp(subAdjGraph);
ldlGolden, schurCGolden = Laplacians.approxChol(llmat, 1)
@test Laplacians.isapprox(ldlSub, ldlGolden);
@test isapprox(schurCSub, schurCGolden);
#also test the condition number, should be 1
cn = abs(Laplacians.condNumber(subAdjGraph, ldlSub, schurCSub,verbose=true));
@test abs(cn -1)<0.01

#test approx chol in a subset but with more ports
llmat = Laplacians.LLmatp(adjGraph);
ldl, schurC = Laplacians.approxChol(llmat, 8);
subset = [1, 3, 5, 7, 9];
ldlSub, schurCSub = Laplacians.subsetApproxChol(ldl, schurC, subset);
#should have the same results with single copy
subAdjGraph = adjGraph[subset, subset];
llmat = Laplacians.LLmatp(subAdjGraph);
ldlGolden, schurCGolden = Laplacians.approxChol(llmat, 4)
@test Laplacians.isapprox(ldlSub, ldlGolden);
@test isapprox(schurCSub, schurCGolden);
#also test the condition number, should be 1
cn = abs(Laplacians.condNumber(subAdjGraph, ldlSub, schurCSub,verbose=true));
@test abs(cn -1)<0.01

#test a graph with 0 deg in ports
#modify icm10post, create node 6
II = [1, 2, 1, 5, 2, 3, 2, 4, 3, 4, 4, 5];
JJ = [2, 1, 5, 1, 3, 2, 4, 2, 4, 3, 5, 4];
VV = [1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2];
adjGraph = SparseArrays.sparse(II, JJ, VV, 6, 6);
#run approximate cholesky
llmat = Laplacians.LLmatp(adjGraph);
ldl, schurC = Laplacians.approxChol(llmat, 3);
#get the results in given subset
subset = [1, 2, 3, 4, 5];
ldlSub, schurCSub = Laplacians.subsetApproxChol(ldl, schurC, subset);
#should have the same results with single copy
subAdjGraph = adjGraph[subset, subset];
llmat = Laplacians.LLmatp(subAdjGraph);
ldlGolden, schurCGolden = Laplacians.approxChol(llmat, 2)
@test Laplacians.isapprox(ldlSub, ldlGolden);
@test isapprox(schurCSub, schurCGolden);
#also test the condition number, should be 1
cn = abs(Laplacians.condNumber(subAdjGraph, ldlSub, schurCSub,verbose=true));
@test abs(cn -1)<0.01

#test a graph with 0 deg inside ports
II = [1, 2, 1, 6, 2, 3, 2, 4, 3, 4, 4, 6];
JJ = [2, 1, 6, 1, 3, 2, 4, 2, 4, 3, 6, 4];
VV = [1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2];
adjGraph = SparseArrays.sparse(II, JJ, VV, 6, 6);
#run approximate cholesky
llmat = Laplacians.LLmatp(adjGraph);
ldl, schurC = Laplacians.approxChol(llmat, 3);
#get the results in given subset
subset = [1, 2, 3, 4, 6];
ldlSub, schurCSub = Laplacians.subsetApproxChol(ldl, schurC, subset);
#should have the same results with single copy
subAdjGraph = adjGraph[subset, subset];
llmat = Laplacians.LLmatp(subAdjGraph);
ldlGolden, schurCGolden = Laplacians.approxChol(llmat, 2)
@test Laplacians.isapprox(ldlSub, ldlGolden);
@test isapprox(schurCSub, schurCGolden);
#also test the condition number, should be 1
cn = abs(Laplacians.condNumber(subAdjGraph, ldlSub, schurCSub,verbose=true));
@test abs(cn -1)<0.01


#now test partitioned functions 
#testcase 1, simple H tree partitioned into two parts
#part 1:
# global indexing                         local indexing
#   1                                       1
#   |                                       |
#   |                                       |
#   3---6                                   3---4
#   |                                       |
#   |                                       |
#   2                                       2
I1 = [1, 3, 2, 3, 3, 4]; 
J1 = [3, 1, 3, 2, 4, 3];
V1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
adjGraph1 = SparseArrays.sparse(I1, J1, V1, 4, 4);

#part 2
#   4                                      1
#   |                                      |
#   |                                      |
#   6                                      3
#   |                                      |
#   |                                      |
#   5                                      2
I2 = [1, 3, 2, 3]; 
J2 = [3, 1, 3, 2];
V2 = [1.0, 1.0, 1.0, 1.0];
adjGraph2 = SparseArrays.sparse(I2, J2, V2, 3, 3);

las = [Laplacians.lap(adjGraph1), Laplacians.lap(adjGraph2)];
portVec1 = [1];
portVec2 = [1];
portVecs = [[1], [1]];
indexOffsets = [0, 3];
numInternalNode = 5;
numInternalNodes = [3, 2];
numPort = 1;
y = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
Laplacians.LaplacianVectorMult!(las, y, numInternalNodes, indexOffsets, portVecs, numInternalNode, numPort);
#y should be all 0
@test sum(abs.(y))<1e-10


#the aggregated graph
#  1  4
#  |  |
#  |  |
#  3--6
#  |  |
#  |  |
#  2  5
II = [1, 3, 2, 3, 4, 6, 5, 6, 3, 6];
JJ = [3, 1, 3, 2, 6, 4, 6, 5, 6, 3];
VV = [1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1];
adjGraph = SparseArrays.sparse(II, JJ, VV, 6, 6);

la = Laplacians.lap(adjGraph);
for ii in 1:10
    x = rand(6);
    y = la*x;
    Laplacians.LaplacianVectorMult!(las, x, numInternalNodes, indexOffsets, portVecs, numInternalNode, numPort);
    #println(isapprox(x,y));
    @test isapprox(x,y);
end

#test the schur complement aggregation
schurC1 = SparseArrays.sparse([], [], Float64[], 1, 1);
schurC2 = SparseArrays.sparse([], [], Float64[], 1, 1);
schurCs = [schurC1, schurC2];
schurC = Laplacians.schurComplement(schurCs, portVecs, numPort);

#now test the condition number
llmat1 = Laplacians.LLmatp(adjGraph1);
llmat2 = Laplacians.LLmatp(adjGraph2);
ldl1, schurC1 = Laplacians.approxChol(llmat1, 1);
ldl2, schurC2 = Laplacians.approxChol(llmat2, 1);
ldls = [ldl1, ldl2];
schurCs = [schurC1, schurC2];
adjGraphs = [adjGraph1, adjGraph2];
@test isapprox(abs(Laplacians.condNumber(adjGraphs, ldls, schurCs, portVecs, numPort, verbose=true)),1);

#let's also test whether aggregated adjGraph is correct
adjGraphT = Laplacians.fullAdjGraph(adjGraphs, portVecs, numPort);
@test isapprox(LinearAlgebra.norm(adjGraph - adjGraphT), 0);

#test how early stop works together with distributed approx chol
llmat1 = Laplacians.LLmatp(adjGraph1);
#add one port for part 1
ldlP1, schurCP1, P1 = Laplacians.approxChol(llmat1, 1, 1);
pAdjGraph1 = adjGraph1[P1, P1];
llmat2 = Laplacians.LLmatp(adjGraph2);
#add one port for part 2
ldlP2, schurCP2, P2 = Laplacians.approxChol(llmat2, 1, 1);
pAdjGraph2 = adjGraph2[P2, P2];
pAdjGraphs = [pAdjGraph1, pAdjGraph2];
ldlPs = [ldlP1, ldlP2];
schurCPs = [schurCP1, schurCP2];
pPortVecs = [[1, 3], [2, 3]]
pNumPort = 3;
cn = abs(Laplacians.condNumber(pAdjGraphs, ldlPs, schurCPs, pPortVecs, pNumPort, verbose = true));
@test isapprox(cn, 1);

#Second testcase grid2_3
#part 1
#global view                local view
# 6----7                    3----4
# |    |                    |    |
# |    |                    |    |
# 1----5                    1----2
I1 = [1, 2, 1, 3, 3, 4, 2, 4];
J1 = [2, 1, 3, 1, 4, 3, 4, 2];
V1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
adjGraph1 = SparseArrays.sparse(I1, J1, V1, 4, 4);
portVec1 = [1, 2, 3];
llmat1 = Laplacians.LLmatp(adjGraph1);
ldl1, schurC1 = Laplacians.approxChol(llmat1, 3);

#part 2
#global view                local view
# 7----8                    3----4
#      |                         |
#      |                         |
# 5----2                    2----1
I2 = [1, 2, 1, 4, 3, 4];
J2 = [2, 1, 4, 1, 4, 3];
V2 = [1.0, 1, 1, 1, 1, 1];
adjGraph2 = SparseArrays.sparse(I2, J2, V2, 4, 4);
portVec2 = [1, 3, 4];
llmat2 = Laplacians.LLmatp(adjGraph2);
ldl2, schurC2 = Laplacians.approxChol(llmat2, 3);


#part 3
#global view                local view
# 3----9                    1----4
# |    |                    |    |
# |    |                    |    |
# 6    7                    2    3
I3 = [1, 2, 1, 4, 3, 4];
J3 = [2, 1, 4, 1, 4, 3];
V3 = [1.0, 1, 1, 1, 1, 1];
adjGraph3 = SparseArrays.sparse(I3, J3, V3, 4, 4);
portVec3 = [2, 3, 5];
llmat3 = Laplacians.LLmatp(adjGraph3);
ldl3, schurC3 = Laplacians.approxChol(llmat3, 3);


#part 4
#global view                local view
# 9----4                    3----1
#      |                         |
#      |                         |
#      8                         2
I4 = [1, 3, 2, 1];
J4 = [3, 1, 1, 2];
V4 = [1.0, 1.0, 1.0, 1.0];
adjGraph4 = SparseArrays.sparse(I4, J4, V4, 3, 3);
portVec4 = [4, 5];
llmat4 = Laplacians.LLmatp(adjGraph4);
ldl4, schurC4 = Laplacians.approxChol(llmat4, 2);

#first calculate required dims and parameters
adjGraphs = [adjGraph1, adjGraph2, adjGraph3, adjGraph4];
las = [Laplacians.lap(adjGraph1), Laplacians.lap(adjGraph2), Laplacians.lap(adjGraph3), Laplacians.lap(adjGraph4)];
portVecs = [portVec1, portVec2, portVec3, portVec4];
numPart = size(adjGraphs, 1);
numNodes = Array{Int64}(undef, numPart);
numPorts = Array{Int64}(undef, numPart);
for ii in 1:numPart
    numNodes[ii] = adjGraphs[ii].m;#adjGraphs[ii].m should be equal to adjGraphs[ii].n
    numPorts[ii] = size(portVecs[ii], 1);       
end
numInternalNodes = numNodes .- numPorts;
indexOffsets = accumulate(+, [0; numInternalNodes[1:end-1]]);
numInternalNode = sum(numInternalNodes);
numPort = 5;

#test laplacian matrix vector multiplication
adjGraph = Laplacians.grid2(3);
(II, JJ, VV) = SparseArrays.findnz(adjGraph);
#permute
P = [1, 5, 2, 6, 7, 8, 3, 9, 4];
IIp = P[II];
JJp = P[JJ];
adjGraphP = SparseArrays.sparse(IIp, JJp, VV, 9, 9);
la = Laplacians.lap(adjGraphP);
for ii in 1:10
    x = rand(9);
    y = la*x;
    Laplacians.LaplacianVectorMult!(las, x, numInternalNodes, indexOffsets, portVecs, numInternalNode, numPort);
    #println(isapprox(x,y));
    @test isapprox(x,y);
end
#let's also test whether aggregated adjGraph is correct
adjGraphT = Laplacians.fullAdjGraph(adjGraphs, portVecs, numPort);
@test isapprox(LinearAlgebra.norm(adjGraphP - adjGraphT), 0);

ldls = [ldl1, ldl2, ldl3, ldl4];
schurCs = [schurC1, schurC2, schurC3, schurC4];
#then test the aggregated schur complement
schurC = Laplacians.schurComplement(schurCs, portVecs, numPort);

@test isapprox(abs(Laplacians.condNumber(adjGraphs, ldls, schurCs, portVecs, numPort, verbose=true)),1);

#Let's test grid_2(3) in distributed processing with early stop
#add one port for part 1
llmat1 = Laplacians.LLmatp(adjGraph1);
ldlP1, schurCP1, P1 = Laplacians.approxChol(llmat1, 3, 1);
pAdjGraph1 = adjGraph1[P1, P1];
pPortVec1 = vcat(1, portVec1.+4);
#add one port for part 2
llmat2 = Laplacians.LLmatp(adjGraph2);
ldlP2, schurCP2, P2 = Laplacians.approxChol(llmat2, 3, 1);
pAdjGraph2 = adjGraph2[P2, P2];
pPortVec2 = vcat(2, portVec2.+4);
#add one port for part 3
llmat3 = Laplacians.LLmatp(adjGraph3);
ldlP3, schurCP3, P3 = Laplacians.approxChol(llmat3, 3, 1);
pAdjGraph3 = adjGraph3[P3, P3];
pPortVec3 = vcat(3, portVec3.+4);
#add one port for part 4
llmat4 = Laplacians.LLmatp(adjGraph4);
ldlP4, schurCP4, P4 = Laplacians.approxChol(llmat4, 2, 1);
pAdjGraph4 = adjGraph4[P4, P4];
pPortVec4 = vcat(4, portVec4.+4);

pAdjGraphs = [pAdjGraph1, pAdjGraph2, pAdjGraph3, pAdjGraph4];
ldlPs = [ldlP1, ldlP2, ldlP3, ldlP4];
schurCPs = [schurCP1, schurCP2, schurCP3, schurCP4];
pPortVecs = [pPortVec1, pPortVec2, pPortVec3, pPortVec4];
pNumPort = 9;
cn = abs(Laplacians.condNumber(pAdjGraphs, ldlPs, schurCPs, pPortVecs, pNumPort, verbose = true));
@test isapprox(cn, 1);


#let's test permuted icm10 post
#part 1
# global view                    local view
# 1---5                          1---4
# |                              |
# |                              |
# 2---4                          2---3
I1 = [1, 2, 1, 4, 2, 3];
J1 = [2, 1, 4, 1, 3, 2];
V1 = [1.0, 1, 1, 1, 2, 2];
adjGraph1 = SparseArrays.sparse(I1, J1, V1, 4, 4);
portVec1 = [1, 2];

#part 2
# global view                    local view
# 5---3                          3---1
# |   |                          |   |
# |   |                          |   |
# ----4                          ----2
I2 = [1, 2, 1, 3, 2, 3];
J2 = [2, 1, 3, 1, 3, 2];
V2 = [1.0, 1, 1, 1, 1, 1];
adjGraph2 = SparseArrays.sparse(I2, J2, V2, 3, 3);
#[1,2] or [2, 1] should be both okay
portVec2 = [2, 1];

adjGraphs = [adjGraph1, adjGraph2];
portVecs = [portVec1, portVec2];

llmat1 = Laplacians.LLmatp(adjGraph1);
llmat2 = Laplacians.LLmatp(adjGraph2);
ldl1, schurC1 = Laplacians.approxChol(llmat1, 2);
ldl2, schurC2 = Laplacians.approxChol(llmat2, 2);
ldls = [ldl1, ldl2];
schurCs = [schurC1, schurC2];
numPort = 2;
@test isapprox(abs(Laplacians.condNumber(adjGraphs, ldls, schurCs, portVecs, numPort, verbose=true)),1);

#lets test larger grid 2's 

nGrid = 1;
#generate one part of grid2(nGrid);
adjGraph1 = Laplacians.grid2(nGrid);
(I1, J1, V1) = SparseArrays.findnz(adjGraph1);
#add edges to bottom line
for ii in 1:nGrid
    append!(I1, ii);
    append!(J1, nGrid*nGrid + 1 + ii);
    append!(V1, 1);

    append!(J1, ii);
    append!(I1, nGrid*nGrid + 1 + ii);
    append!(V1, 1);
end
#add edges to vertical line
for ii in 1:nGrid
    append!(I1, (ii-1)*nGrid + 1);
    append!(J1, nGrid*nGrid + 1 + nGrid + ii);
    append!(V1, 1);

    append!(J1, (ii-1)*nGrid + 1);
    append!(I1, nGrid*nGrid + 1 + nGrid + ii);
    append!(V1, 1);
end
#add the bottom line
for ii in 1:nGrid
    append!(I1, nGrid*nGrid + ii);
    append!(J1, nGrid*nGrid + ii + 1);
    append!(V1, 1);

    append!(J1, nGrid*nGrid + ii);
    append!(I1, nGrid*nGrid + ii + 1);
    append!(V1, 1);
end
adjGraphi = SparseArrays.sparse(I1, J1, V1, (nGrid+1)*(nGrid+1), (nGrid+1)*(nGrid+1));
portVec1 = collect(1:(2*nGrid +1));
portVec2 = [1;collect((nGrid+2): (3*nGrid+1))];
portVec3 = [1;collect((2nGrid+2): (4*nGrid+1))];
portVec4 = [1;collect((3nGrid+2): (4*nGrid+1)); collect(2: (nGrid+1))];

adjGraphs = [adjGraphi, adjGraphi, adjGraphi, adjGraphi];
las = [Laplacians.lap(adjGraphi), Laplacians.lap(adjGraphi), Laplacians.lap(adjGraphi), Laplacians.lap(adjGraphi)];
portVecs = [portVec1, portVec2, portVec3, portVec4];
numPart = size(adjGraphs, 1);
numNodes = Array{Int64}(undef, numPart);
numPorts = Array{Int64}(undef, numPart);
for ii in 1:numPart
    numNodes[ii] = adjGraphs[ii].m;#adjGraphs[ii].m should be equal to adjGraphs[ii].n
    numPorts[ii] = size(portVecs[ii], 1);       
end
numInternalNodes = numNodes .- numPorts;
indexOffsets = accumulate(+, [0; numInternalNodes[1:end-1]]);
numInternalNode = sum(numInternalNodes);
numPort = 5;


#test the Laplacian matrix vector multiplication
adjGraph = Laplacians.grid2(3);
(II, JJ, VV) = SparseArrays.findnz(adjGraph);
#permute
P = [3, 9, 4, 8, 5, 6, 2, 7, 1];
IIp = P[II];
JJp = P[JJ];
adjGraphP = SparseArrays.sparse(IIp, JJp, VV, 9, 9);
la = Laplacians.lap(adjGraphP);
for ii in 1:10
    x = rand(9);
    y = la*x;
    Laplacians.LaplacianVectorMult!(las, x, numInternalNodes, indexOffsets, portVecs, numInternalNode, numPort);
    #println(isapprox(x,y));
    @test isapprox(x,y);
end

#let's also test whether aggregated adjGraph is correct
adjGraphT = Laplacians.fullAdjGraph(adjGraphs, portVecs, numPort);
@test isapprox(LinearAlgebra.norm(adjGraphP - adjGraphT), 0);

#test condition number
llmati = Laplacians.LLmatp(adjGraphi);
ldli, schurCi = Laplacians.approxChol(llmati, 3);

ldls = [ldli, ldli, ldli, ldli];
schurCs = [schurCi, schurCi, schurCi, schurCi];
adjGraphs = [adjGraphi, adjGraphi, adjGraphi, adjGraphi];
@test isapprox(abs(Laplacians.condNumber(adjGraphs, ldls, schurCs, portVecs, numPort, verbose=true)),1);

#test larger grid2
nGrid = 50;
#generate one part of grid2(nGrid);
adjGraph1 = Laplacians.grid2(nGrid);
(I1, J1, V1) = SparseArrays.findnz(adjGraph1);
#add edges to bottom line
for ii in 1:nGrid
    append!(I1, ii);
    append!(J1, nGrid*nGrid + 1 + ii);
    append!(V1, 1);

    append!(J1, ii);
    append!(I1, nGrid*nGrid + 1 + ii);
    append!(V1, 1);
end
#add edges to vertical line
for ii in 1:nGrid
    append!(I1, (ii-1)*nGrid + 1);
    append!(J1, nGrid*nGrid + 1 + nGrid + ii);
    append!(V1, 1);

    append!(J1, (ii-1)*nGrid + 1);
    append!(I1, nGrid*nGrid + 1 + nGrid + ii);
    append!(V1, 1);
end
#add the bottom line
for ii in 1:nGrid
    append!(I1, nGrid*nGrid + ii);
    append!(J1, nGrid*nGrid + ii + 1);
    append!(V1, 1);

    append!(J1, nGrid*nGrid + ii);
    append!(I1, nGrid*nGrid + ii + 1);
    append!(V1, 1);
end
adjGraphi = SparseArrays.sparse(I1, J1, V1, (nGrid+1)*(nGrid+1), (nGrid+1)*(nGrid+1));
portVec1 = collect(1:(2*nGrid +1));
portVec2 = [1;collect((nGrid+2): (3*nGrid+1))];
portVec3 = [1;collect((2nGrid+2): (4*nGrid+1))];
portVec4 = [1;collect((3nGrid+2): (4*nGrid+1)); collect(2: (nGrid+1))];

adjGraphs = [adjGraphi, adjGraphi, adjGraphi, adjGraphi];
las = [Laplacians.lap(adjGraphi), Laplacians.lap(adjGraphi), Laplacians.lap(adjGraphi), Laplacians.lap(adjGraphi)];
portVecs = [portVec1, portVec2, portVec3, portVec4];
numPart = size(adjGraphs, 1);
numNodes = Array{Int64}(undef, numPart);
numPorts = Array{Int64}(undef, numPart);
for ii in 1:numPart
    numNodes[ii] = adjGraphs[ii].m;#adjGraphs[ii].m should be equal to adjGraphs[ii].n
    numPorts[ii] = size(portVecs[ii], 1);       
end
numInternalNodes = numNodes .- numPorts;
indexOffsets = accumulate(+, [0; numInternalNodes[1:end-1]]);
numInternalNode = sum(numInternalNodes);
numPort = 4*nGrid + 1;

#test condition number
llmati = Laplacians.LLmatp(adjGraphi);
ldli, schurCi = Laplacians.approxChol(llmati, 2*nGrid+1);

ldls = [ldli, ldli, ldli, ldli];
schurCs = [schurCi, schurCi, schurCi, schurCi];
@test abs(Laplacians.condNumber(adjGraphs, ldls, schurCs, portVecs, numPort, verbose=true))<20;#should be 12-15

#even larger
nGrid = 1500;
#generate one part of grid2(nGrid);
adjGraph1 = Laplacians.grid2(nGrid);
(I1, J1, V1) = SparseArrays.findnz(adjGraph1);
#add edges to bottom line
for ii in 1:nGrid
    append!(I1, ii);
    append!(J1, nGrid*nGrid + 1 + ii);
    append!(V1, 1);

    append!(J1, ii);
    append!(I1, nGrid*nGrid + 1 + ii);
    append!(V1, 1);
end
#add edges to vertical line
for ii in 1:nGrid
    append!(I1, (ii-1)*nGrid + 1);
    append!(J1, nGrid*nGrid + 1 + nGrid + ii);
    append!(V1, 1);

    append!(J1, (ii-1)*nGrid + 1);
    append!(I1, nGrid*nGrid + 1 + nGrid + ii);
    append!(V1, 1);
end
#add the bottom line
for ii in 1:nGrid
    append!(I1, nGrid*nGrid + ii);
    append!(J1, nGrid*nGrid + ii + 1);
    append!(V1, 1);

    append!(J1, nGrid*nGrid + ii);
    append!(I1, nGrid*nGrid + ii + 1);
    append!(V1, 1);
end
adjGraphi = SparseArrays.sparse(I1, J1, V1, (nGrid+1)*(nGrid+1), (nGrid+1)*(nGrid+1));
portVec1 = collect(1:(2*nGrid +1));
portVec2 = [1;collect((nGrid+2): (3*nGrid+1))];
portVec3 = [1;collect((2nGrid+2): (4*nGrid+1))];
portVec4 = [1;collect((3nGrid+2): (4*nGrid+1)); collect(2: (nGrid+1))];

adjGraphs = [adjGraphi, adjGraphi, adjGraphi, adjGraphi];
las = [Laplacians.lap(adjGraphi), Laplacians.lap(adjGraphi), Laplacians.lap(adjGraphi), Laplacians.lap(adjGraphi)];
portVecs = [portVec1, portVec2, portVec3, portVec4];
numPart = size(adjGraphs, 1);
numNodes = Array{Int64}(undef, numPart);
numPorts = Array{Int64}(undef, numPart);
for ii in 1:numPart
    numNodes[ii] = adjGraphs[ii].m;#adjGraphs[ii].m should be equal to adjGraphs[ii].n
    numPorts[ii] = size(portVecs[ii], 1);       
end
numInternalNodes = numNodes .- numPorts;
indexOffsets = accumulate(+, [0; numInternalNodes[1:end-1]]);
numInternalNode = sum(numInternalNodes);
numPort = 4*nGrid + 1;

#test condition number
llmati = Laplacians.LLmatp(adjGraphi);
ldli, schurCi = Laplacians.approxChol(llmati, 2*nGrid+1);

ldls = [ldli, ldli, ldli, ldli];
schurCs = [schurCi, schurCi, schurCi, schurCi];
@test abs(Laplacians.condNumber(adjGraphs, ldls, schurCs, portVecs, numPort, verbose=true))<40;#should be 20-25

#test dump LDL, schurC and load LDL, schurC
adjGraph = Laplacians.grid2(3);
llmat = Laplacians.LLmatp(adjGraph);
ldl, schurC = Laplacians.approxChol(llmat, 1);
#dump the obtained ldl, schurC
Laplacians.dumpApproxFact("grid_2_3_1", ldl, schurC);
#read the dumped file
ldlBak, schurCBak = Laplacians.loadApproxFact("grid_2_3_1");
@test Laplacians.isapprox(ldl, ldlBak);
@test isapprox(schurC, schurCBak);
rm("grid_2_3_1")

#test dump and load with larger testcases
adjGraph = Laplacians.grid2(100);
llmat = Laplacians.LLmatp(adjGraph);
ldl, schurC = Laplacians.approxChol(llmat, 50);
#dump the obtained ldl, schurC
Laplacians.dumpApproxFact("grid_2_100_50", ldl, schurC);
#read the dumped file
ldlBak, schurCBak = Laplacians.loadApproxFact("grid_2_100_50");
@test Laplacians.isapprox(ldl, ldlBak);
@test isapprox(schurC, schurCBak);
rm("grid_2_100_50")

println("End of testByExport")
