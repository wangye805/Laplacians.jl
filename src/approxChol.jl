#=

approxChol Laplacian solver by Daniel A. Spielman, 2017.
This algorithm is an implementation of an approximate edge-by-edge elimination
algorithm inspired by the Approximate Gaussian Elimination algorithm of
Kyng and Sachdeva.

For usage exaples, see http://danspielman.github.io/Laplacians.jl/latest/usingSolvers/index.html

There are two versions of this solver:
one that fixes the order of elimination beforehand,
and one that adapts the order to eliminate verties of low degree.
These use different data structures.
LLOrdMat is for the fixed order, and LLmatp is for the adaptive order.

These coes produce a structure we call LDLinv that is then used in the solve.
The structure of this code is as follows:

The data structures appear in approxCholTypes.jl
We then have the outline:

* constructors for LLmatp and LLMatOrd
* get_ll_col and compress_ll_col : used inside the elimination
* approxChol : the main routine
* LDLsolver, and its forward and backward solve the apply LDLinv
* approxchol_lap: the main solver, which calls approxchol_lap1 on connected
    components.
    This then calls one of approxchol_lapWdeg, approxchol_lapGiven or approxchol_lapGreedy,
    depending on the parameters.

* approxchol_lapChol - for producing a Cholesky factor instead of an LDLinv.
  might be useful if optimized.
* data structures that are used for the adaptive low-degree version to
  choose the next vertex.

=#

"""
    params = ApproxCholParams(order, output)
order can be one of
* :deg (by degree, adaptive),
* :wdeg (by original wted degree, nonadaptive),
* :given
"""
mutable struct ApproxCholParams
    order::Symbol
    stag_test::Integer
end

ApproxCholParams() = ApproxCholParams(:deg, 5)
ApproxCholParams(sym::Symbol) = ApproxCholParams(sym, 5)

LDLinv(a::SparseMatrixCSC{Tval,Tind}) where {Tind,Tval} =
  LDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n),Tind[],Tval[],zeros(Tval,a.n))

LDLinv(a::LLMatOrd{Tind,Tval}) where {Tind,Tval} =
  LDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n),Tind[],Tval[],zeros(Tval,a.n))

LDLinv(a::LLmatp{Tind,Tval}) where {Tind,Tval} =
  LDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n),Tind[],Tval[],zeros(Tval,a.n))

LDL(a::LLmatp{Tind,Tval}, np::Tind) where {Tind,Tval} = LDL(a.n-np, np, 1, Tind[], 1, Tind[], Tind[], Tval[], 1, Tind[], Tind[], Tval[], Tval[]);

function LLmatp(a::SparseMatrixCSC{Tval,Tind}) where {Tind,Tval}
    n = size(a,1)
    m = nnz(a)

    degs = zeros(Tind,n)

    flips = flipIndex(a)

    cols = Array{LLp{Tind,Tval}}(undef, n)
    llelems = Array{LLp{Tind,Tval}}(undef, m)

    @inbounds for i in 1:n
        degs[i] = a.colptr[i+1] - a.colptr[i]
        if degs[i]==0
            continue;
        end
        ind = a.colptr[i]
        j = a.rowval[ind]
        v = a.nzval[ind]
        llpend = LLp{Tind,Tval}(j,v)
        next = llelems[ind] = llpend
        for ind in (a.colptr[i]+one(Tind)):(a.colptr[i+1]-one(Tind))
            j = a.rowval[ind]
            v = a.nzval[ind]
            next = llelems[ind] = LLp{Tind,Tval}(j,v,next)
        end
        cols[i] = next
    end

    @inbounds for i in 1:n
        if degs[i]==0
            continue;
        end
        for ind in a.colptr[i]:(a.colptr[i+1]-one(Tind))
            llelems[ind].reverse = llelems[flips[ind]]
        end
    end

    return LLmatp{Tind,Tval}(n, degs, cols, llelems)
end

"""
Print all column in an LLmatp matrix
"""
function print_ll_mat(llmat::LLmatp)
    @inbounds for i in 1:llmat.n
        println("column $(i)");
        print_ll_col(llmat, i)
    end
end

"""
  Print a column in an LLmatp matrix.
  This is here for diagnostics.
"""
function print_ll_col(llmat::LLmatp, i::Int)
    ll = llmat.cols[i]
    println("col $i, row $(ll.row) : $(ll.val)")

    while ll.next != ll
        ll = ll.next
        println("col $i, row $(ll.row) : $(ll.val)")
    end
end

function LLMatOrd(a::SparseMatrixCSC{Tval,Tind}) where {Tind,Tval}
    n = size(a,1)
    m = nnz(a)

    cols = zeros(Tind, n)
    llelems = Array{LLord{Tind,Tval}}(undef, m)

    ptr = one(Tind)

    @inbounds for i in Tind(1):Tind(n-1)
        next = zero(Tind)

        for ind in (a.colptr[i]):(a.colptr[i+1]-one(Tind))
            j = a.rowval[ind]
            if (i < j)

              v = a.nzval[ind]
              llelems[ptr] = LLord{Tind,Tval}(j, next, v)
              next = ptr
              ptr += one(Tind)

            end
        end
        cols[i] = next
    end

    return LLMatOrd{Tind,Tval}(n, cols, llelems)
end

function LLMatOrd(a::SparseMatrixCSC{Tval,Tind}, perm::Array) where {Tind,Tval}
    n = size(a,1)
    m = nnz(a)

    invp = invperm(perm)

    cols = zeros(Tind, n)
    llelems = Array{LLord{Tind,Tval}}(undef, m)

    ptr = one(Tind)

    @inbounds for i0 in Tind(1):Tind(n)
        i = invp[i0]
        next = zero(Tind)

        for ind in (a.colptr[i0]):(a.colptr[i0+1]-one(Tind))
            j = invp[a.rowval[ind]]
            if (i < j)

              v = a.nzval[ind]
              llelems[ptr] = LLord{Tind,Tval}(j, next, v)
              next = ptr
              ptr += one(ptr)

            end
        end
        cols[i] = next
    end

    return LLMatOrd{Tind,Tval}(n, cols, llelems)
end

"""
  Print a column in an LLMatOrd matrix.
  This is here for diagnostics.
"""
function print_ll_col(llmat::LLMatOrd, i::Int)
    ptr = llmat.cols[i]
    while ptr != 0
      ll = llmat.lles[ptr]
      println("col $i, row $(ll.row) : $(ll.val)")

      ptr = ll.next
    end
end



#=============================================================

The approximate factorization

=============================================================#

function get_ll_col(llmat::LLmatp{Tind,Tval},
  i,
  colspace::Vector{LLp{Tind,Tval}}) where {Tind,Tval}


    ll = llmat.cols[i]
    len = 0
    @inbounds while ll.next != ll

        if ll.val > zero(Tval)
            len = len+1
            if (len > length(colspace))
                push!(colspace,ll)
            else
                colspace[len] = ll
            end
        end

        ll = ll.next
    end

    if ll.val > zero(Tval)
        len = len+1
        if (len > length(colspace))
            push!(colspace,ll)
        else
            colspace[len] = ll
        end
    end

    return len
end

function get_ll_col(llmat::LLMatOrd{Tind,Tval},
  i,
  colspace::Vector{LLcol{Tind,Tval}}) where {Tind,Tval}

    ptr = llmat.cols[i]
    len = 0
    @inbounds while ptr != 0

        #if ll.val > 0
            len = len+1

            # should not be an lles - is an abuse
            item = LLcol(llmat.lles[ptr].row, ptr, llmat.lles[ptr].val)
            if (len > length(colspace))
                push!(colspace,item)
            else
                colspace[len] = item
            end
        #end

        ptr = llmat.lles[ptr].next
    end

    return len
end
#specially used for schur complement
function compressCol!(a::LLmatp{Tind,Tval},
  colspace::Vector{LLp{Tind,Tval}},
  len::Int;
  forSchurC=true) where {Tind,Tval}

    o = Base.Order.ord((x, y)->(x[1]<y[1]||(x[1]==y[1]&&x[2]<y[2])), x->(x.row, x.val), false, Base.Order.Forward)

    sort!(colspace, 1, len, QuickSort, o)

    ptr = 0
    currow::Tind = 0

    c = colspace

    @inbounds for i in 1:len

        if c[i].row != currow
            currow = c[i].row
            ptr = ptr+1
            c[ptr] = c[i]

        else
            c[ptr].val = c[ptr].val + c[i].val
            if (!forSchurC) #for schur complement, we do not need the next line
                c[i].reverse.val = zero(Tval)
            end
        end
    end


    o = Base.Order.ord(isless, x->x.val, false, Base.Order.Forward)
    sort!(colspace, 1, ptr, QuickSort, o)

    return ptr
end


function compressCol!(a::LLmatp{Tind,Tval},
  colspace::Vector{LLp{Tind,Tval}},
  len::Int,
  pq::ApproxCholPQ{Tind}) where {Tind,Tval}

    o = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)

    sort!(colspace, 1, len, QuickSort, o)

    ptr = 0
    currow::Tind = 0

    c = colspace

    @inbounds for i in 1:len

        if c[i].row != currow
            currow = c[i].row
            ptr = ptr+1
            c[ptr] = c[i]

        else
            c[ptr].val = c[ptr].val + c[i].val
            c[i].reverse.val = zero(Tval)

            approxCholPQDec!(pq, currow)
        end
    end


    o = Base.Order.ord(isless, x->x.val, false, Base.Order.Forward)
    sort!(colspace, 1, ptr, QuickSort, o)

    return ptr
end

function compressCol!(
  colspace::Vector{LLcol{Tind,Tval}},
  len::Int
  ) where {Tind,Tval}

    o = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)

    sort!(colspace, one(len), len, QuickSort, o)

    c = colspace

    ptr = 0
    currow = c[1].row
    curval = c[1].val
    curptr = c[1].ptr

    @inbounds for i in 2:len

        if c[i].row != currow

            ptr = ptr+1
            c[ptr] = LLcol(currow, curptr, curval)  # next is abuse here: reall keep where it came from.

            currow = c[i].row
            curval = c[i].val
            curptr = c[i].ptr

        else

            curval = curval + c[i].val

        end

    end

    # emit the last row

    ptr = ptr+1
    c[ptr] = LLcol(currow, curptr, curval)

    o = Base.Order.ord(isless, x->x.val, false, Base.Order.Forward)
    sort!(colspace, one(ptr), ptr, QuickSort, o)

    return ptr
end


function approxChol(a::LLMatOrd{Tind,Tval}) where {Tind,Tval}
    n = a.n

    # need to make custom one without col info later.
    ldli = LDLinv(a)
    ldli_row_ptr = one(Tind)

    d = zeros(Tval,n)

    colspace = Array{LLcol{Tind,Tval}}(undef, n)
    cumspace = Array{Tval}(undef, n)
    #vals = Array(Tval,n) # will be able to delete this

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)


    for i in Tind(1):Tind(n-1)

        ldli.col[i] = i  # will get rid of this with new data type
        ldli.colptr[i] = ldli_row_ptr

        len = get_ll_col(a, i, colspace)

        len = compressCol!(colspace, len)

        csum = zero(Tval)
        for ii in 1:len
            #vals[ii] = colspace[ii].val    # if immut, no need for vals
            csum = csum + colspace[ii].val
            cumspace[ii] = csum
        end
        wdeg = csum

        colScale = one(Tval)

        for joffset in 1:(len-1)

            llcol = colspace[joffset]
            w = llcol.val * colScale
            j = llcol.row

            f = w/(wdeg)

            #vals[joffset] = zero(Tval)

            # kind = Laplacians.blockSample(vals,k=1)[1]
            r = rand() * (csum - cumspace[joffset]) + cumspace[joffset]
            koff = searchsortedfirst(cumspace,r,one(len),len,o)

            k = colspace[koff].row

            newEdgeVal = w*(one(Tval)-f)

            # create edge (j,k) with newEdgeVal
            # do it by reassigning ll
            if j < k # put it in col j
                jhead = a.cols[j]
                a.lles[llcol.ptr] = LLord(k, jhead, newEdgeVal)
                #ll.next = jhead
                #ll.val = newEdgeVal
                #ll.row = k
                a.cols[j] = llcol.ptr
            else # put it in col k
              khead = a.cols[k]
              a.lles[llcol.ptr] = LLord(j, khead, newEdgeVal)
              #ll.next = khead
              #ll.val = newEdgeVal
              #ll.row = j
              a.cols[k] = llcol.ptr
            end

            colScale = colScale*(one(Tval)-f)
            #wdeg = wdeg*(1.0-f)^2
            wdeg = wdeg - 2w + w^2/wdeg

            push!(ldli.rowval,j)
            push!(ldli.fval, f)
            ldli_row_ptr = ldli_row_ptr + one(Tind)

            # push!(ops, IJop(i,j,1-f,f))  # another time suck


        end # for joffset


        llcol = colspace[len]
        w = llcol.val * colScale
        j = llcol.row

        push!(ldli.rowval,j)
        push!(ldli.fval, one(Tval))
        ldli_row_ptr = ldli_row_ptr + one(Tind)

        d[i] = w

    end # for i


    ldli.colptr[n] = ldli_row_ptr

    ldli.d = d

    return ldli
end

# this one is greedy on the degree - also a big win
function approxChol(a::LLmatp{Tind,Tval}) where {Tind,Tval}
    n = a.n

    ldli = LDLinv(a)
    ldli_row_ptr = one(Tind)

    d = zeros(n)

    pq = ApproxCholPQ(a.degs)

    it = 1

    colspace = Array{LLp{Tind,Tval}}(undef, n)
    cumspace = Array{Tval}(undef, n)
    vals = Array{Tval}(undef, n) # will be able to delete this

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)

    @inbounds while it < n

        i = approxCholPQPop!(pq)

        ldli.col[it] = i # conversion!
        ldli.colptr[it] = ldli_row_ptr

        it = it + 1

        len = get_ll_col(a, i, colspace)

        len = compressCol!(a, colspace, len, pq)  #3hog

        csum = zero(Tval)
        for ii in 1:len
            vals[ii] = colspace[ii].val
            csum = csum + colspace[ii].val
            cumspace[ii] = csum
        end
        wdeg = csum

        colScale = one(Tval)

        for joffset in 1:(len-1)

            ll = colspace[joffset]
            w = vals[joffset] * colScale
            j = ll.row
            revj = ll.reverse

            f = w/(wdeg)

            vals[joffset] = zero(Tval)

            # kind = Laplacians.blockSample(vals,k=1)[1]
            r = rand() * (csum - cumspace[joffset]) + cumspace[joffset]
            koff = searchsortedfirst(cumspace,r,one(len),len,o)

            k = colspace[koff].row

            approxCholPQInc!(pq, k)

            newEdgeVal = f*(one(Tval)-f)*wdeg

            # fix row k in col j
            revj.row = k   # dense time hog: presumably becaus of cache
            revj.val = newEdgeVal
            revj.reverse = ll

            # fix row j in col k
            khead = a.cols[k]
            a.cols[k] = ll
            ll.next = khead
            ll.reverse = revj
            ll.val = newEdgeVal
            ll.row = j


            colScale = colScale*(one(Tval)-f)
            wdeg = wdeg*(one(Tval)-f)^2

            push!(ldli.rowval,j)
            push!(ldli.fval, f)
            ldli_row_ptr = ldli_row_ptr + one(Tind)

            # push!(ops, IJop(i,j,1-f,f))  # another time suck


        end # for


        ll = colspace[len]
        w = vals[len] * colScale
        j = ll.row
        revj = ll.reverse

        if it < n
            approxCholPQDec!(pq, j)
        end

        revj.val = zero(Tval)

        push!(ldli.rowval,j)
        push!(ldli.fval, one(Tval))
        ldli_row_ptr = ldli_row_ptr + one(Tind)

        d[i] = w

    end

    ldli.colptr[it] = ldli_row_ptr

    ldli.d = d

    return ldli
end

#the function to factorize one node 
#debug:         print the reduce degs and PRNGs used for each node 
#fixedRandom:   use the given PRNGs for each node
#When: 
#1. debug == false, fixedRandom == false, order, PRNGs, reduceDegs not used and empty 
#2. debug == true, fixedRandom == false, order not used and empty, PRNGs, reduceDegs will bereturned as an output tuple
#3. debug == false, fixedRandom == true, order, PRNGs provided as input, reduceDegs empty and not used
#4. debug == true, fixedRandom == true, order, PRNGs provided as input, reduceDegs returned as output variable
function factorizeOneNode(a::LLmatp{Tind,Tval},
                          i::Tind, 
                          colspace::Array{LLp{Tind,Tval}}, 
                          cumspace::Array{Tval},
                          vals::Array{Tval},
                          len::Tind,
                          csum::Tval,
                          ldl::LDL{Tind, Tval},
                          pq::ApproxCholPQ{Tind},
                          debug::Bool, 
                          fixedRandom::Bool, 
                          reduceDegs::Array{Tind, 1}, 
                          PRNGs::Array{Array{Tval,1}, 1}) where {Tind,Tval}
    addColumn!(ldl, i);

    #add the diagonal element
    addDiag!(ldl, csum)

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)
    colScale = one(Tval)
    wdeg = csum;   
    #record reduce degs for each iteration of reduction
    if(debug)
        reduceDegs[i] = len;
        if(!fixedRandom)
            PRNGs[i] = Tval[];
        end
    end
    #first len-1 edges are special
    for joffset in 1:(len-1)

        ll = colspace[joffset]
        w = vals[joffset] * colScale
        j = ll.row
        revj = ll.reverse

        #println("deal with node $(j) in col $(i)");

        f = w/(wdeg)

        addOffDiag!(ldl, j, vals[joffset])

        vals[joffset] = zero(Tval)
        if(fixedRandom)
            randVal = PRNGs[i][joffset];
        else
            randVal = rand();
            if debug
                push!(PRNGs[i], randVal);
            end
        end
        # kind = Laplacians.blockSample(vals,k=1)[1]
        r = randVal * (csum - cumspace[joffset]) + cumspace[joffset]
        koff = searchsortedfirst(cumspace,r,one(len),len,o)

        k = colspace[koff].row
        if(!fixedRandom)
            approxCholPQInc!(pq, k)
        end
        newEdgeVal = f*(one(Tval)-f)*wdeg
        #println("newEdgeVal is $(newEdgeVal)");

        # fix row k in col j
        revj.row = k   # dense time hog: presumably becaus of cache
        revj.val = newEdgeVal
        revj.reverse = ll

        # fix row j in col k
        khead = a.cols[k]
        a.cols[k] = ll
        ll.next = khead
        ll.reverse = revj
        ll.val = newEdgeVal
        ll.row = j


        colScale = colScale*(one(Tval)-f)
        wdeg = wdeg*(one(Tval)-f)^2
        
    end # for

    ll = colspace[len]
    w = vals[len] * colScale
    j = ll.row
    revj = ll.reverse

    #println("deal with last node $(j) in col $(i)");
    revj.val = zero(Tval)

    addOffDiag!(ldl,j, vals[len]);

    if(!fixedRandom)
        approxCholPQDec!(pq, j)
    end


end

#the approximate cholesy with early stop
#nReduce is the actual node to be reduced
function approxChol(a::LLmatp{Tind,Tval}, nPorts::Tind, numAddPorts::Tind) where {Tind,Tval}
    #dim of the original graph
    nGraph = a.n
    #num of nodes pushed into priority queue
    n = a.n - nPorts
    
    #initialize an empty LDL
    ldl = LDL(a, nPorts)
    
    #now only feed in first n degs
    pq = ApproxCholPQ(a.degs, n)

    colspace = Array{LLp{Tind,Tval}}(undef, nGraph)
    cumspace = Array{Tval}(undef, nGraph)
    vals = Array{Tval}(undef, nGraph) # will be able to delete this

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)

    it = 1
    nReduce = nGraph - nPorts - numAddPorts;
    @inbounds while it <= nReduce
        i = approxCholPQPop!(pq)
        addColumn!(ldl, i);

        it = it + 1

        len = get_ll_col(a, i, colspace)

        len = compressCol!(a, colspace, len, pq)  #3hog

        csum = zero(Tval)
        for ii in 1:len
            vals[ii] = colspace[ii].val
            csum = csum + colspace[ii].val
            cumspace[ii] = csum
        end
        wdeg = csum

        #add the diagonal element
        addDiag!(ldl, csum)

        colScale = one(Tval)
        
        for joffset in 1:(len-1)

            ll = colspace[joffset]
            w = vals[joffset] * colScale
            j = ll.row
            revj = ll.reverse

            f = w/(wdeg)

            addOffDiag!(ldl, j, vals[joffset])

            vals[joffset] = zero(Tval)
            randVal = rand();

            r = randVal * (csum - cumspace[joffset]) + cumspace[joffset]
            koff = searchsortedfirst(cumspace,r,one(len),len,o)

            k = colspace[koff].row
            #println("k is $(k)");
            approxCholPQInc!(pq, k)

            newEdgeVal = f*(one(Tval)-f)*wdeg
            #println("newEdgeVal is $(newEdgeVal)");

            # fix row k in col j
            revj.row = k   # dense time hog: presumably becaus of cache
            revj.val = newEdgeVal
            revj.reverse = ll

            # fix row j in col k
            khead = a.cols[k]
            a.cols[k] = ll
            ll.next = khead
            ll.reverse = revj
            ll.val = newEdgeVal
            ll.row = j

            colScale = colScale*(one(Tval)-f)
            wdeg = wdeg*(one(Tval)-f)^2
            
        end # for

        ll = colspace[len]
        w = vals[len] * colScale
        j = ll.row
        revj = ll.reverse

        #println("deal with last node $(j) in col $(i)");

        addOffDiag!(ldl,j, vals[len]);
        revj.val = zero(Tval)
    end
    #pop all unReduced nodes as ports
    while it<=n
        i = approxCholPQPop!(pq);
        addPort!(ldl, i);
        it+=1;
    end
    finishColumn!(ldl);
    reform!(ldl);
    P = [ldl.col;n+1:a.n];
    schurC = schurComplement!(a, nReduce, numAddPorts, ldl.col);
    #make ldl.col 1:ldl.n
    ldl.col = collect(1:ldl.n);
    return ldl, schurC, P
end


#the approximate cholesky with parameter of number of ports (>=1)
#debug:         print the reduce degs and PRNGs used for each node 
#fixedRandom:   use the given PRNGs for each node
#When: 
#1. debug == false, fixedRandom == false, order, PRNGs, reduceDegs not used and empty 
#2. debug == true, fixedRandom == false, order not used and empty, PRNGs, reduceDegs will bereturned as an output tuple
#3. debug == false, fixedRandom == true, order, PRNGs provided as input, reduceDegs empty and not used
#4. debug == true, fixedRandom == true, order, PRNGs provided as input, reduceDegs returned as output variable
function approxChol(a::LLmatp{Tind,Tval}, nPorts::Tind; debug = false, fixedRandom = false, order::Array{Tind, 1} = Tind[], PRNGs::Array{Array{Tval,1}, 1} = Array{Tval, 1}[]) where {Tind,Tval}
    nGraph = a.n
    n = a.n - nPorts
    
    #replace with our LDL class
    ldl = LDL(a, nPorts)

    #now only feed in first n degs
    #not acutally used when fixed random is true
    pq = ApproxCholPQ(a.degs, n)

    it = 1

    colspace = Array{LLp{Tind,Tval}}(undef, nGraph)
    cumspace = Array{Tval}(undef, nGraph)
    vals = Array{Tval}(undef, nGraph) # will be able to delete this

    if(debug)
        #record the reduced degs at each step
        reduceDegs = Array{Int64}(undef, n);
        if(!fixedRandom)
            PRNGs = Array{Array{Tval,1}, 1}(undef, n);
        end
    else
        reduceDegs = Tind[];
        if(!fixedRandom)
            PRNGs = Array{Tval, 1}[];
        end
    end

    @inbounds while it <= n
        if(fixedRandom)
            #use given order
            i = order[it];
        else
            i = approxCholPQPop!(pq)
        end
        #println("pop out node $(i)");

        it = it + 1;

        #pre-processing edge weights for node i
        len = get_ll_col(a, i, colspace);
        if(!fixedRandom)
            len = compressCol!(a, colspace, len, pq)  #3hog
        else
            len = compressCol!(a, colspace, len; forSchurC=false)  #3hog
        end

        csum = zero(Tval)
        for ii in 1:len
            vals[ii] = colspace[ii].val
            csum = csum + colspace[ii].val
            cumspace[ii] = csum
        end

        #factorize node i
        factorizeOneNode(a, i, colspace, cumspace, vals, len, csum, ldl, pq,
                         debug, 
                         fixedRandom, 
                         reduceDegs, 
                         PRNGs);

    end
    finishColumn!(ldl)  
    #change the ldl to be the canonical form
    canonicalForm!(ldl)
    #post processing for debugging info
    debugInfo = [];
    if debug
        debugInfo = [reduceDegs, PRNGs];
    end
    #also the remaining schur complement graph (as adjacency matrix)
    return ldl, schurComplement!(a, n), debugInfo
end

#the function to dump LDL
#for debugging only
#index has already been unified with c/c++ format to start from 0
#can choose not to dump value by verbose option
function dumpLDL(fh::IO, ldl::LDL{Tind, Tval}, verbose = true) where {Tind,Tval}
    #print n
    write(fh, "n: \n$(ldl.n)\n")
    #print m
    write(fh, "m: \n$(ldl.m)\n")
    #print perm_idx
    write(fh, "perm_idx: \n$(ldl.perm_idx - 1)\n")
    #print col array
    write(fh, "col: \n")
    for idx in ldl.col
        write(fh, "$(idx-1)\n")
    end
    #print current_row_ptr_A
    write(fh, "current_row_ptr_A: \n$(ldl.current_row_ptr_A - 1)\n");
    #print colptr_A
    write(fh, "colptr_A: \n")
    for idx in ldl.colptr_A
        write(fh, "$(idx-1)\n")
    end
    #print rowval_A
    write(fh, "rowval_A: \n")
    for idx in ldl.rowval_A
        write(fh, "$(idx - 1)\n")
    end
    if(verbose)
        #print nzval_A
        write(fh, "nzval_A: \n")
        for val in ldl.nzval_A
            write(fh, "$(val)\n")
        end
    end

    #print current_row_ptr_B
    write(fh, "current_row_ptr_B: \n$(ldl.current_row_ptr_B - 1)\n");
    #print colptr_B
    write(fh, "colptr_B: \n")
    for idx in ldl.colptr_B
        write(fh, "$(idx-1)\n")
    end
    #print rowval_B
    write(fh, "rowval_B: \n")
    for idx in ldl.rowval_B
        write(fh, "$(idx - 1)\n")
    end
    if(verbose)
        #print nzval_B
        write(fh, "nzval_B: \n")
        for val in ldl.nzval_B
            write(fh, "$(val)\n")
        end
    end
    #print diagA if necessary
    if(verbose)
        write(fh, "diag: \n")
        for val in ldl.diagA
            write(fh, "$(val)\n")
        end
    end
end

#the function to compare two LDL's
function isapprox(ldl1::LDL{Tind, Tval}, ldl2::LDL{Tind, Tval}) where {Tind,Tval}
    #compare each field of LDL
    #compare m
    if(!Base.isapprox(ldl1.n, ldl2.n))
        return false;
    end
    #compare n
    if(!Base.isapprox(ldl1.m, ldl2.m))
        return false;
    end
    #compare col
    if(!Base.isapprox(ldl1.col, ldl2.col))
        return false;
    end
    #compare perm_idx
    if(!Base.isapprox(ldl1.perm_idx, ldl2.perm_idx))
        return false;
    end
    #compare current_row_ptr_A
    if(!Base.isapprox(ldl1.current_row_ptr_A, ldl2.current_row_ptr_A))
        return false;
    end
    #compare colptr_A
    if(!Base.isapprox(ldl1.colptr_A, ldl2.colptr_A))
        return false;
    end
    #compare rowval_A
    if(!Base.isapprox(ldl1.rowval_A, ldl2.rowval_A))
        return false;
    end
    #compare nzval_A
    if(!Base.isapprox(ldl1.nzval_A, ldl2.nzval_A))
        return false;
    end
    #compare current_row_ptr_B
    if(!Base.isapprox(ldl1.current_row_ptr_B, ldl2.current_row_ptr_B))
        return false;
    end
    #compare colptr_B
    if(!Base.isapprox(ldl1.colptr_B, ldl2.colptr_B))
        return false;
    end
    #compare rowval_B
    if(!Base.isapprox(ldl1.rowval_B, ldl2.rowval_B))
        return false;
    end
    #compare nzval_B
    if(!Base.isapprox(ldl1.nzval_B, ldl2.nzval_B))
        return false;
    end
    #compare diagA
    if(!Base.isapprox(ldl1.diagA, ldl2.diagA))
        return false;
    end
    return true;
end

#the function to dump SchurComplement
#for debugging only
#index has already been unified with c/c++ format to start from 0
#can choose not to dump value by verbose option
function dumpSchurC(fh::IO, schurC::SparseMatrixCSC{Tval, Tind}, verbose = true) where {Tind,Tval}
    #print m
    write(fh, "m: \n$(schurC.m)\n");
    #print n
    write(fh, "n: \n$(schurC.n)\n");
    #print colptr
    write(fh, "colptr: \n")
    for idx in schurC.colptr
        write(fh, "$(idx-1)\n")
    end
    #print nnz
    write(fh, "nnz: \n$(size(schurC.rowval, 1))\n");
    #print rowval
    write(fh, "rowval: \n")
    for idx in schurC.rowval
        write(fh, "$(idx-1)\n")
    end
    if(verbose)
        #print nzval
        write(fh, "nzval: \n")
        for val in schurC.nzval
            write(fh, "$(val)\n")
        end
    end
end

#the function to dump approximate factorization results
#for debugging only
#index has already been unified with C/C++ to start from 0
#can choose to dump connection only or also with value by verbose
function dumpApproxFact(filename::AbstractString, ldl::LDL{Tind, Tval}, schurC::SparseMatrixCSC{Tval, Tind}, verbose = true) where {Tind,Tval}
    fh = open(filename,"w")
    #print LDL first
    dumpLDL(fh, ldl, verbose);
    #then print schurC
    dumpSchurC(fh, schurC, verbose);
    close(fh)
end

# fix the pseudo-randomness, and allow partial factorization
# nPorts denote the number of ports, >=1 (because GND is always a port node)
# PRNG are given random numbers, usually of size m (# of edges) for the worst case
function approxChol(a::LLmatp{Tind,Tval}, nPorts::Tind, PRNG::Array{Tval}) where {Tind,Tval}
    exactFactorization=true
    nGraph = a.n
    n = a.n - nPorts
    
    #replace with our LDL class
    ldl = LDL(a, nPorts)

    #now only feed in first n degs
    pq = ApproxCholPQ(a.degs, n)

    it = 1

    colspace = Array{LLp{Tind,Tval}}(undef, nGraph)
    cumspace = Array{Tval}(undef, nGraph)
    vals = Array{Tval}(undef, nGraph) # will be able to delete this

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)

    idxPRNG = 1;

    @inbounds while it <= n

        i = approxCholPQPop!(pq)
        #println("pop out node $(i)");

        #PrintApproxCholPQ(pq);

        addColumn!(ldl, i);

        it = it + 1

        len = get_ll_col(a, i, colspace)

        len = compressCol!(a, colspace, len, pq)  #3hog

        csum = zero(Tval)
        for ii in 1:len
            vals[ii] = colspace[ii].val
            csum = csum + colspace[ii].val
            cumspace[ii] = csum
        end
        wdeg = csum

        #add the diagonal element
        addDiag!(ldl, csum)

        colScale = one(Tval)
        
        if len>=3
            exactFactorization = false;
        end

        for joffset in 1:(len-1)

            ll = colspace[joffset]
            w = vals[joffset] * colScale
            j = ll.row
            revj = ll.reverse

            #println("deal with node $(j) in col $(i)");

            f = w/(wdeg)

            addOffDiag!(ldl, j, vals[joffset])

            vals[joffset] = zero(Tval)

            # kind = Laplacians.blockSample(vals,k=1)[1]
            #r = rand() * (csum - cumspace[joffset]) + cumspace[joffset]
            r = PRNG[idxPRNG] * (csum - cumspace[joffset]) + cumspace[joffset]
            idxPRNG = idxPRNG + 1;
            koff = searchsortedfirst(cumspace,r,one(len),len,o)

            k = colspace[koff].row
            #println("k is $(k)");
            approxCholPQInc!(pq, k)

            newEdgeVal = f*(one(Tval)-f)*wdeg
            #println("newEdgeVal is $(newEdgeVal)");

            # fix row k in col j
            revj.row = k   # dense time hog: presumably becaus of cache
            revj.val = newEdgeVal
            revj.reverse = ll

            # fix row j in col k
            khead = a.cols[k]
            a.cols[k] = ll
            ll.next = khead
            ll.reverse = revj
            ll.val = newEdgeVal
            ll.row = j


            colScale = colScale*(one(Tval)-f)
            wdeg = wdeg*(one(Tval)-f)^2
            
        end # for

        ll = colspace[len]
        w = vals[len] * colScale
        j = ll.row
        revj = ll.reverse

        #println("deal with last node $(j) in col $(i)");

        addOffDiag!(ldl,j, vals[len]);

        if it < nGraph
            approxCholPQDec!(pq, j)
        end

        revj.val = zero(Tval)
        #lets print the graph in llmat 
        #print_ll_mat(a);
    end

    finishColumn!(ldl)  
    #println("is exact factorization? $(exactFactorization)")
    #also the remaining schur complement graph (as adjacency matrix)
    return ldl, schurComplement!(a, n)
end

#=============================================================
The function to get the remaining schur complement from linked list matrix
n is the number of internal nodes
=============================================================#
function schurComplement!(a::LLmatp{Tind,Tval}, n::Tind) where {Tind, Tval}
    I = Tind[];
    J = Tind[];
    V = Tval[];
    #prepare for loading each column
    colspace = Array{LLp{Tind,Tval}}(undef, a.n)

    @inbounds for ii in n+1:1:a.n
        if a.degs[ii]==0
            continue
        end
        len = get_ll_col(a, ii, colspace)
        len = compressCol!(a, colspace, len)  #3hog
        @inbounds for jj in 1:1:len
            #println("colspace[$(jj)]: row-> $(colspace[jj].row), val->$(colspace[jj].val)");
            if colspace[jj].row>n
                #(i, j, v)
                push!(I, colspace[jj].row - n);
                push!(J, ii - n);
                push!(V, colspace[jj].val);
            end
        end
    end
    return sparse(I, J, V, a.n-n, a.n-n);
end

#=============================================================
The function to get the remaining schur complement from linked list matrix from the early stop mode
=============================================================#
function schurComplement!(a::LLmatp{Tind,Tval}, n::Tind, numAddPorts::Tind, permVec::Array{Tind, 1}) where {Tind, Tval}
    #construct the perm mapping from permVec
    permMap = Array{Tind, 1}(undef, a.n);
    for ii=1:1:a.n
        if(ii<=size(permVec,1))
            permMap[permVec[ii]] = ii;
        else
            permMap[ii] = ii;
        end
    end
    I = Tind[];
    J = Tind[];
    V = Tval[];
    #prepare for loading each column
    colspace = Array{LLp{Tind,Tval}}(undef, a.n)

    @inbounds for ii in n+1:1:a.n
        if a.degs[ii]==0
            continue
        end
        if(ii<=size(permVec,1))
            len = get_ll_col(a, permVec[ii], colspace)
        else
            len = get_ll_col(a, ii, colspace)
        end
        len = compressCol!(a, colspace, len)  #3hog
        @inbounds for jj in 1:1:len
            push!(I, permMap[colspace[jj].row] - n);
            push!(J, ii - n);
            push!(V, colspace[jj].val);
        end
    end
    return sparse(I, J, V, a.n-n, a.n-n);
end



#=============================================================
The methods for new LDL
=============================================================#
function addColumn!(ldl::LDL{Tind, Tval}, ii::Tind)where {Tind,Tval}
    push!(ldl.col, ii);

    push!(ldl.colptr_A, ldl.current_row_ptr_A);
    push!(ldl.colptr_B, ldl.current_row_ptr_B);
end

function addDiag!(ldl::LDL{Tind, Tval}, csum::Tval) where {Tind,Tval}
    push!(ldl.diagA, csum);
    ldl.perm_idx += 1;
end

function addOffDiag!(ldl::LDL{Tind, Tval}, jj::Tind, val::Tval) where {Tind,Tval}
    if(jj<=ldl.n)#add to part A
        push!(ldl.rowval_A, jj);
        push!(ldl.nzval_A, val);
        ldl.current_row_ptr_A += 1;
    else#add to partB
        push!(ldl.rowval_B, jj-ldl.n);
        push!(ldl.nzval_B, val);
        ldl.current_row_ptr_B += 1;
    end
end

function finishColumn!(ldl::LDL{Tind, Tval}) where {Tind,Tval}
    push!(ldl.colptr_A, ldl.current_row_ptr_A);
    push!(ldl.colptr_B, ldl.current_row_ptr_B);
    #convert raw row index in LA part into permuted index
    #ldl.col is always of size ldl.n even if in early stop
    permMap = Array{Tind, 1}(undef, ldl.n)
    @inbounds for ii in 1:ldl.n
        permMap[ldl.col[ii]] = ii;
    end
    @inbounds for ii in 1:(length(ldl.rowval_A))
        ldl.rowval_A[ii] = permMap[ldl.rowval_A[ii]];
    end
end

#re-format LA, LB part as the canonical form of sparse Matrix CSC
function canonicalForm!(ldl::LDL{Tind, Tval}) where {Tind,Tval}
    #deal with LA part first
    @inbounds for ii=1:1:(size(ldl.colptr_A,1)-1)
        rowvals = ldl.rowval_A[ldl.colptr_A[ii]:(ldl.colptr_A[ii+1]-1)];
        nzvals = ldl.nzval_A[ldl.colptr_A[ii]:(ldl.colptr_A[ii+1]-1)];
        p = sortperm(rowvals);
        #write the sorted rowval and nzval
        ldl.rowval_A[ldl.colptr_A[ii]:(ldl.colptr_A[ii+1]-1)] = rowvals[p];
        ldl.nzval_A[ldl.colptr_A[ii]:(ldl.colptr_A[ii+1]-1)] = nzvals[p];
    end
    #then deal with LB part
    @inbounds for ii=1:1:(size(ldl.colptr_B,1)-1)
        rowvals = ldl.rowval_B[ldl.colptr_B[ii]:(ldl.colptr_B[ii+1]-1)];
        nzvals = ldl.nzval_B[ldl.colptr_B[ii]:(ldl.colptr_B[ii+1]-1)];
        p = sortperm(rowvals);
        #write the sorted rowval and nzval
        ldl.rowval_B[ldl.colptr_B[ii]:(ldl.colptr_B[ii+1]-1)] = rowvals[p];
        ldl.nzval_B[ldl.colptr_B[ii]:(ldl.colptr_B[ii+1]-1)] = nzvals[p];
    end
end

#add additional ports into LDL
function addPort!(ldl::LDL{Tind, Tval}, ii::Tind) where {Tind,Tval}
    push!(ldl.col, ii);
    ldl.perm_idx += 1; 
end

#change the ldl dimension due to early stop
function reform!(ldl::LDL{Tind, Tval}) where {Tind,Tval}
    numAddPorts = size(ldl.col, 1) - size(ldl.diagA, 1);
    n = ldl.n - numAddPorts;
    m = ldl.m + numAddPorts;
    #transform the ldl into the canonical form so that we can have valid sparse matrix csc
    canonicalForm!(ldl);
    LA = SparseArrays.SparseMatrixCSC{Tval, Tind}(ldl.n, n, ldl.colptr_A, ldl.rowval_A, ldl.nzval_A);
    LB = SparseArrays.SparseMatrixCSC{Tval, Tind}(ldl.m, n, ldl.colptr_B, ldl.rowval_B, ldl.nzval_B);
    L = [LA;LB];
    #re-distribute LA and LB according to new dim
    LA = L[1:n, :];
    LB = L[n+1:end, :];
    #update data members
    ldl.n = n;
    ldl.m = m;
    ldl.current_row_ptr_A = size(LA.rowval,1) + 1;
    ldl.colptr_A = LA.colptr;
    ldl.rowval_A = LA.rowval;
    ldl.nzval_A = LA.nzval;

    ldl.current_row_ptr_B = size(LB.rowval,1) + 1;
    ldl.colptr_B = LB.colptr;
    ldl.rowval_B = LB.rowval;
    ldl.nzval_B = LB.nzval;
end

#=============================================================

The routines that do the solve.

=============================================================#

function LDLsolver(ldli::LDLinv, b::Vector)
    y = copy(b)

    forward!(ldli, y)

    @inbounds for i in 1:(length(ldli.d))
        if ldli.d[i] != 0
            y[i] /= ldli.d[i]
        end
    end

    backward!(ldli, y)

    mu = mean(y)
    @inbounds for i in eachindex(y)
        y[i] = y[i] - mu
    end

    return y
end


function forward!(ldli::LDLinv{Tind,Tval}, y::Vector) where {Tind,Tval}

    @inbounds for ii in 1:length(ldli.col)
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-one(Tind)

        yi = y[i]

        for jj in j0:(j1-1)
            j = ldli.rowval[jj]
            y[j] += ldli.fval[jj] * yi
            yi *= (one(Tval)-ldli.fval[jj])
        end
        j = ldli.rowval[j1]
        y[j] += yi
        y[i] = yi
    end
end

function backward!(ldli::LDLinv{Tind,Tval}, y::Vector) where {Tind,Tval}
    o = one(Tind)
    @inbounds for ii in length(ldli.col):-1:1
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-o

        j = ldli.rowval[j1]
        yi = y[i]
        yi = yi + y[j]

        for jj in (j1-o):-o:j0
            j = ldli.rowval[jj]
            yi = (one(Tval)-ldli.fval[jj])*yi + ldli.fval[jj]*y[j]
        end
        y[i] = yi
    end
end

#forwardSolve API required by Intel Pardiso package
#calculate b2-B^TA^{-1}b1
function forwardSolve!(ldl::LDL{Tind, Tval}, y::Vector)where {Tind, Tval}
    b1 = y[1:ldl.n];
    b2 = y[(ldl.n+1):(ldl.n+ldl.m)];
    
    #now we want to get b2-LB*LA^-1P^Tb1
    #P^Tb_1
    permute!(b1, ldl.col);
    #LA^-1
    #b1 = LA\b1;
    invNLA!(ldl, b1);
    #b2-LB*y
    #b2 = b2 - LB*b1;
    b2 = b2-nLB(ldl, b1);
    #re-assign b2 portion
    y[(ldl.n+1):(ldl.n+ldl.m)] = b2;
end

#normalized LA solve: nLA x = y
function invNLA!(ldl::LDL{Tind, Tval}, y::Vector) where {Tind, Tval}
    @inbounds for i in 1:ldl.n-1
        @inbounds for j in ldl.colptr_A[i]:(ldl.colptr_A[i+1]-1)
            y[ldl.rowval_A[j]] = y[ldl.rowval_A[j]] + ldl.nzval_A[j]/ldl.diagA[i]*y[i];
        end
    end
end

#normalized LA' solve nLA'^{-1}y        
function invNLAT!(ldl::LDL{Tind, Tval}, y::Vector) where {Tind, Tval}
    @inbounds for i in ldl.n-1:-1:1
        @inbounds for j in ldl.colptr_A[i]:(ldl.colptr_A[i+1]-1)
            y[i] = y[i] + ldl.nzval_A[j]/ldl.diagA[i]*y[ldl.rowval_A[j]];
        end
    end
end

#matrix vector multiplication nLB*x
function nLB(ldl::LDL{Tind, Tval}, x::Vector) where {Tind, Tval}
    y = zeros(Tval, ldl.m);
    @inbounds for i in 1:1:ldl.n
        @inbounds for j in ldl.colptr_B[i]:(ldl.colptr_B[i+1]-1)
            y[ldl.rowval_B[j]] = y[ldl.rowval_B[j]] -ldl.nzval_B[j]/ldl.diagA[i]*x[i];
        end
    end
    return y;
end

#matrix vector multiplication nLB'*x
function nLBT(ldl::LDL{Tind, Tval}, x::Vector) where {Tind, Tval}
    y = zeros(Tval, ldl.n);
    @inbounds for i in 1:1:ldl.n
        @inbounds for j in ldl.colptr_B[i]:(ldl.colptr_B[i+1]-1)
            y[i] = y[i] - ldl.nzval_B[j]/ldl.diagA[i]*x[ldl.rowval_B[j]];
        end
    end
    return y;
end

#the function to do y = invD*x
function invD!(ldl::LDL{Tind, Tval}, x::Vector) where {Tind, Tval}
    @inbounds for i in 1:1:ldl.n
        x[i] = x[i]/ldl.diagA[i];
    end
end

#backwardSolve API required by Intel Pardiso package
#calculate b2-B^TA^{-1}b1
function backwardSolve!(ldl::LDL{Tind, Tval}, y::Vector)where {Tind, Tval}
    b1 = y[1:ldl.n];
    b2 = y[(ldl.n+1):(ldl.n+ldl.m)];#this is x2 part
    
    #we want to get P*LA^{-T}(invD*LA^{-1}*P^T*b1 - LB^Tb2)
    #P^T(b_1)
    permute!(b1, ldl.col);
    #invD*LA^-1
    #b1 = invD*(LA\b1);
    invNLA!(ldl, b1);
    invD!(ldl, b1);

    #LB^T*b2
    #b2 = (LB')*b2; 
    b2 = nLBT(ldl, b2);
    #b1 = (LA')\(b1-b2);
    b1 = b1 - b2;
    invNLAT!(ldl, b1);
    #P(x)
    invpermute!(b1, ldl.col);
    
    #re-assign b1 portion
    y[1:ldl.n] = b1;
end

function solve!(ldl::LDL{Tind, Tval}, y::Vector, s::SparseMatrixCSC{Tval, Tind})where {Tind, Tval}
    forwardSolve!(ldl, y);
    #schur^{-1}
    S = lap(s)[1:end-1, 1:end-1];
    if ldl.m>1
        x2 = [S\y[ldl.n+1:ldl.n+ldl.m-1]; 0];
    else
        x2 = [0];
    end
    y[(ldl.n+1):(ldl.n+ldl.m)] = x2;
    backwardSolve!(ldl,y);
end


#=
  An attempt at an efficient solver for the case when y is a matrix.
  Have not yet found a meaningful speedup

function LDLsolver(ldli::LDLinv, b::Matrix)
    y = copy(b)

    (d, n) = size(y)
    @assert n == length(ldli.col)+1

    forward!(ldli, y)

    @inbounds for i in 1:(length(ldli.d))
        if ldli.d[i] != 0
          @simd for j in 1:d
            y[j,i] = y[j,i] / ldli.d[i]
          end
        end
    end

    backward!(ldli, y)

    @inbounds for j in 1:size(y,1)
        mu = mean(y[j,:])

        for i in 1:size(y,2)
            y[j,i] = y[j,i] - mu
        end
    end

    return y
end



function forward!{Tind,Tval}(ldli::LDLinv{Tind,Tval}, y::Matrix)

    (d, n) = size(y)
    @assert n == length(ldli.col)+1

    #yi = zeros(y[:,1])

    @inbounds for ii in 1:length(ldli.col)
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-one(Tind)

        for jj in j0:(j1-1)
            j = ldli.rowval[jj]
            @simd for k in 1:d
              y[k,j] = y[k,j] + (ldli.fval[jj] * y[k,i])
              y[k,i] = y[k,i] * (one(Tval)-ldli.fval[jj])
          end
        end
        j = ldli.rowval[j1]

        @simd for k in 1:d
          y[k,j] = y[k,j] + y[k,i]
        end
    end
end

function backward!{Tind,Tval}(ldli::LDLinv{Tind,Tval}, y::Matrix)
    o = one(Tind)

    (d, n) = size(y)
    @assert n == length(ldli.col)+1

    yi = zeros(y[:,1])

    @inbounds for ii in length(ldli.col):-1:1
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-o

        j = ldli.rowval[j1]
        #copy!(yi, y[:,i])

        @simd for k in 1:d
          y[k,i] = y[k,i] + y[k,j]
        end

        for jj in (j1-o):-o:j0
            j = ldli.rowval[jj]
            @simd for k in 1:d
              y[k,i] = (one(Tval)-ldli.fval[jj])*y[k,i] + ldli.fval[jj].*y[k,j]
            end
        end
        #y[:,i] = yi
    end
end

=#


"""
    solver = approxchol_lap(a); x = solver(b);
    solver = approxchol_lap(a; tol::Real=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams())

A heuristic by Daniel Spielman inspired by the linear system solver in https://arxiv.org/abs/1605.02353 by Rasmus Kyng and Sushant Sachdeva.  Whereas that paper eliminates vertices one at a time, this eliminates edges one at a time.  It is probably possible to analyze it.
The `ApproxCholParams` let you choose one of three orderings to perform the elimination.

* ApproxCholParams(:given) - in the order given.
    This is the fastest for construction the preconditioner, but the slowest solve.
* ApproxCholParams(:deg) - always eliminate the node of lowest degree.
    This is the slowest build, but the fastest solve.
* ApproxCholParams(:wdeg) - go by a perturbed order of wted degree.

For more info, see http://danspielman.github.io/Laplacians.jl/latest/usingSolvers/index.html
"""
function approxchol_lap(a::SparseMatrixCSC{Tv,Ti};
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams()) where {Tv,Ti}

  if minimum(a.nzval) < 0
      error("Adjacency matrix can not have negative edge weights")
  end

    return Laplacians.lapWrapComponents(approxchol_lap1, a,
    verbose=verbose,
    tol=tol,
    maxits=maxits,
    maxtime=maxtime,
    pcgIts=pcgIts,
    params=params)


end

function approxchol_lapGreedy(a::SparseMatrixCSC;
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams())

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts

  t1 = time()

  la = lap(a) # a hit !?

  llmat = LLmatp(a)
  ldli = approxChol(llmat)
  F(b) = LDLsolver(ldli, b)

  if verbose
    println("Using greedy degree ordering. Factorization time: ", time()-t1)
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
  end

  if verbose
      println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
  end


  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b .- mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)

end

function approxchol_lapGiven(a::SparseMatrixCSC;
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams())

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts

  t1 = time()

  la = lap(a)

  llmat = LLMatOrd(a)
  ldli = approxChol(llmat)
  F(b) = LDLsolver(ldli, b)

  if verbose
    println("Using given ordering. Factorization time: ", time()-t1)
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
  end

  if verbose
      println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
  end


  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b .- mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)

end

function approxchol_lapWdeg(a::SparseMatrixCSC;
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams())

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts

  t1 = time()

  la = lap(a)

  v = vec(sum(a,dims=1))
  v = v .* (1 .+ rand(length(v)))
  p = sortperm(v)

  llmat = LLMatOrd(a,p)
  ldli = approxChol(llmat)

  ip = invperm(p)
  ldlip = LDLinv(p[ldli.col], ldli.colptr, p[ldli.rowval], ldli.fval, ldli.d[ip]);

  F = function(b)
    x = zeros(size(b))
    x = LDLsolver(ldlip, b)
    #x[p] = LDLsolver(ldli, b[p])
    return x
  end

  if verbose
    println("Using wted degree ordering. Factorization time: ", time()-t1)
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
  end

  if verbose
      println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
  end


  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b .- mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)

end



function approxchol_lap1(a::SparseMatrixCSC{Tv,Ti};
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams()) where {Tv,Ti}

    tol_ =tol
    maxits_ =maxits
    maxtime_ =maxtime
    verbose_ =verbose
    pcgIts_ =pcgIts


    if params.order == :deg

      return approxchol_lapGreedy(a,
        verbose=verbose,
        tol=tol,
        maxits=maxits,
        maxtime=maxtime,
        pcgIts=pcgIts,
        params=params)


    elseif params.order == :wdeg

      return approxchol_lapWdeg(a,
        verbose=verbose,
        tol=tol,
        maxits=maxits,
        maxtime=maxtime,
        pcgIts=pcgIts,
        params=params)


    else
      return approxchol_lapGiven(a,
        verbose=verbose,
        tol=tol,
        maxits=maxits,
        maxtime=maxtime,
        pcgIts=pcgIts,
        params=params)


    end

end

"""
    solver = approxchol_sddm(sddm); x = solver(b);
    solver = approxchol_sddm(sddm; tol=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams())

Solves sddm systems by wrapping approxchol_lap.
Not yet optimized directly for sddm.

For more info, see http://danspielman.github.io/Laplacians.jl/latest/usingSolvers/index.html
"""
approxchol_sddm = sddmWrapLap(approxchol_lap)




#===============================

  Checking the condition number

=================================#

"""
    cn = condNumber(a, ldli; verbose=false)

Given an adjacency matrix a and an ldli computed by approxChol,
this computes the condition number.
"""
function condNumber(a, ldli; verbose=false)
  la = lap(a)

  # construct the square operator
  g = function(b)

    y = copy(b)

    #=
    mu = mean(y)
    @inbounds for i in eachindex(y)
        y[i] = y[i] - mu
    end
      =#

    @inbounds for i in 1:(length(ldli.d))
        if ldli.d[i] != 0
            y[i] /= (ldli.d[i])^(1/2)
        else
            y[i] = 0
        end
    end

    backward!(ldli, y)

    y = la * y

    forward!(ldli, y)

    @inbounds for i in 1:(length(ldli.d))
        if ldli.d[i] != 0
            y[i] /= (ldli.d[i])^(1/2)
        else
            y[i] = 0
        end
    end

    #=
    mu = mean(y)
    @inbounds for i in eachindex(y)
        y[i] = y[i] - mu
    end
    =#

    return y
  end

  gOp = SqLinOp(true,1.0,size(a,1),g)
  upper = eigs(gOp;nev=1,which=:LM,tol=1e-2)[1][1]

  g2(b) = upper*b - g(b)
  g2Op = SqLinOp(true,1.0,size(a,1),g2)
  lower = upper - eigs(g2Op;nev=2,which=:LM,tol=1e-2)[1][2]

  if verbose
      println("lower: ", lower, ", upper: ", upper);
  end

  return upper/lower

end

"""
condition number for given LDL and schur complement
"""
function condNumber(a, ldl::LDL{Tind, Tval}, schurC::SparseMatrixCSC{Tval, Tind}; verbose=false)where {Tind, Tval}
    @assert((ldl.n+ldl.m)==size(a,1))
    canonicalForm!(ldl);
    la = lap(a)
    # construct the square operator
    g = function(b)
        y = copy(b);
        #the asym op is ldl^-1*la*b
        y = la*y;
        solve!(ldl, y, schurC);
        return y
    end
    gOp = SqLinOp(false,1.0,size(a,1),g)
    upper = eigs(gOp;nev=1,which=:LM,tol=1e-2)[1][1]
    
    g2(b) = upper*b - g(b)
    g2Op = SqLinOp(false,1.0,size(a,1),g2)
    lower = upper - eigs(g2Op;nev=2,which=:LM,tol=1e-2)[1][2]
    
    if verbose
        println("lower: ", lower, ", upper: ", upper);
    end
    
    return upper/lower
end

"""
get the true factor: [LA, 0
                      LB, I] 
"""
function getL(ldl::LDL{Tind, Tval})where {Tind, Tval}
    #make sure that ldl is already in its canonical form
    canonicalForm!(ldl);
    nLA = SparseArrays.SparseMatrixCSC{Tval, Tind}(ldl.n, ldl.n, ldl.colptr_A, ldl.rowval_A, ldl.nzval_A);
    nLB = SparseArrays.SparseMatrixCSC{Tval, Tind}(ldl.m, ldl.n, ldl.colptr_B, ldl.rowval_B, ldl.nzval_B);
    nL = [nLA; nLB];
    #normalize to get tru nL
    for ii=1:1:nL.n
        for jj=nL.colptr[ii]:1:(nL.colptr[ii+1] - 1)
            nL.nzval[jj] = nL.nzval[jj]/sqrt(ldl.diagA[ii]);
        end
    end
    nL = SparseArrays.sparse(1:ldl.n, 1:ldl.n, sqrt.(ldl.diagA), ldl.m + ldl.n, ldl.n) - nL;
    I = SparseArrays.sparse(1:ldl.m, 1:ldl.m, ones(ldl.m));
    E = SparseArrays.sparse(Tind[], Tind[], Tval[], ldl.n, ldl.m);
    L = hcat(nL, [E;I]);
end

"""
Support from LDL, schurC to original graph adj matrix: max eig of la^-1*L*L^T
The asymmetric version is not stable, please use the symmetric version
"""
function support_asym(a, ldl::LDL{Tind, Tval}, schurC::SparseMatrixCSC{Tval, Tind}; verbose=false)where {Tind, Tval}
    #check problem dim
    @assert((ldl.n+ldl.m)==size(a,1))
    #make ldl to be the canonical form
    #getL function below already has canonicalForm enforced.
    #canonicalForm!(ldl);
    L = getL(ldl);
    lc = Laplacians.lap(schurC);
    #the solver for la
    f = Laplacians.approxchol_lap(a, tol=1e-5);
    #define the linear operation for max eigs
    g = function(b)
        y = copy(b);
        #perform PT
        x1 = y[1:ldl.n];
        permute!(x1, ldl.col);
        y[1:ldl.n] = x1;

        #perform L^T*y
        y = L'*y;

        #perform Schur Complement multiplication
        y[ldl.n+1:end] = lc*y[ldl.n+1:end];

        #perform L*y
        y = L*y;

        #perform P y
        x1 = y[1:ldl.n];
        invpermute!(x1, ldl.col);
        y[1:ldl.n] = x1;

        return f(y);
    end
    op = SqLinOp(false,1.0,size(a,1),g)
    return abs(eigs(op;nev=1,which=:LM,tol=1e-2)[1][1]);
end

"""
Calculate the support from ldl, schurC to adjGraph
In other words, cacluate lambda(ldl^-1*la)
"""
function support(ldl::LDL{Tind, Tval}, schurC::SparseMatrixCSC{Tval, Tind}, a::SparseMatrixCSC{Tval, Tind}; verbose=false)where {Tind, Tval}
    @assert((ldl.n+ldl.m)==size(a,1))
    canonicalForm!(ldl);
    la = lap(a)
    F = LinearAlgebra.cholesky(Laplacians.lap(schurC)[1:end-1, 1:end-1]);
    #construct the square operator 
    g = function(b)
        y = copy(b);
        if ldl.m >1
            x2 = [F.UP\y[ldl.n+1:ldl.n+ldl.m-1]; 0];
        else
            x2 = [0];
        end
        y[(ldl.n+1):(ldl.n+ldl.m)] = x2;

        backwardSolve!(ldl, y);
        y = la*y;
        forwardSolve!(ldl, y);

        if ldl.m >1
            x2 = [F.PtL\y[ldl.n+1:ldl.n+ldl.m-1]; 0];
        else
            x2 = [0];
        end
        y[(ldl.n+1):(ldl.n+ldl.m)] = x2;

        return y
    end
    gOp = SqLinOp(true,1.0,size(a,1),g)
    upper = eigs(gOp;nev=1,which=:LM,tol=1e-2)[1][1]
    return upper; 
end

"""
Calculate the support function using the symmetric operand
"""
function support(a, ldl::LDL{Tind, Tval}, schurC::SparseMatrixCSC{Tval, Tind}; verbose=false)where {Tind, Tval}
    #check problem dim
    @assert((ldl.n+ldl.m)==size(a,1))
    #canonical form will be enforced by getL
    #canonicalForm!(ldl);
    L = getL(ldl);
    
    #get the factorization form of lap(schurC)
    F = LinearAlgebra.cholesky(Laplacians.lap(schurC)[1:end-1, 1:end-1]);
    Lc = SparseArrays.sparse(F.L)[invperm(F.p),:];
    LAug = [Lc; -sum(Lc, dims=1)];
    #Laplacian solver
    f = Laplacians.approxchol_lap(a, tol=1e-5);
    #define the symmetric linear function
    g = function(b)
        y = copy(b);

        #perform Lc*y
        x2 = y[ldl.n+1:end-1];
        y[ldl.n+1:end] = LAug*x2;

        #perform L*y
        y = L*y;

        #perform P*y
        x1 = y[1:ldl.n];
        invpermute!(x1, ldl.col);
        y[1:ldl.n] = x1;

        #perform A^-1
        y = f(y);

        #perform P^T y
        x1 = y[1:ldl.n];
        permute!(x1, ldl.col);
        y[1:ldl.n] = x1;

        #perform L^T*y
        y = L'*y;

        #perform Lc^T*y
        x2 = y[ldl.n+1:end];
        y[ldl.n+1:end] = [LAug'*x2; 0];
        
        return y;
    end
    op = SqLinOp(true,1.0,size(a,1),g)
    return eigs(op;nev=1,which=:LM,tol=1e-2)[1][1];
end

"""
condition number for given LDL and schur complement in the symmetric form
symmetric form will be more stable
"""
function condNumberSym(a, ldl::LDL{Tind, Tval}, schurC::SparseMatrixCSC{Tval, Tind}; verbose=false)where {Tind, Tval}
    @assert((ldl.n+ldl.m)==size(a,1))
    canonicalForm!(ldl);
    la = lap(a)
    F = LinearAlgebra.cholesky(Laplacians.lap(schurC)[1:end-1, 1:end-1]);
    #construct the square operator 
    g = function(b)
        y = copy(b);
        if ldl.m >1
            x2 = [F.UP\y[ldl.n+1:ldl.n+ldl.m-1]; 0];
        else
            x2 = [0];
        end
        y[(ldl.n+1):(ldl.n+ldl.m)] = x2;

        backwardSolve!(ldl, y);
        y = la*y;
        forwardSolve!(ldl, y);

        if ldl.m >1
            x2 = [F.PtL\y[ldl.n+1:ldl.n+ldl.m-1]; 0];
        else
            x2 = [0];
        end
        y[(ldl.n+1):(ldl.n+ldl.m)] = x2;

        return y
    end
    gOp = SqLinOp(true,1.0,size(a,1),g)
    upper = eigs(gOp;nev=1,which=:LM,tol=1e-2)[1][1]
    
    g2(b) = upper*b - g(b)
    g2Op = SqLinOp(true,1.0,size(a,1),g2)
    lower = upper - eigs(g2Op;nev=2,which=:LM,tol=1e-2)[1][2]
    
    if verbose
        println("lower: ", lower, ", upper: ", upper);
    end
    
    return upper/lower
end

#obtain the subset of an approximate factorization result
function subsetApproxChol(ldl::LDL{Tind, Tval}, schurC::SparseMatrixCSC{Tval, Tind}, subset::Array{Tind, 1}) where {Tind, Tval}
    #first sort the subset indices
    sort!(subset);
    #create mapping from global index to index in subset
    subsetMap = Dict{Tind, Tind}();
    for ii in 1:1:length(subset)
        subsetMap[subset[ii]]=ii;
    end

    #get to know how many internal/external nodes in subset
    #which is also the dimension of subset approx chol
    numInternalNode = searchsortedlast(subset, ldl.n);
    numPort = length(subset) - numInternalNode;

    #transform the ldl into the canonical form so that we can have valid sparse matrix csc
    canonicalForm!(ldl);
    #create the square matrix [LA, 0; LB, 0]
    LA = SparseArrays.SparseMatrixCSC{Tval, Tind}(ldl.n, ldl.n, ldl.colptr_A, ldl.rowval_A, ldl.nzval_A);
    LB = SparseArrays.SparseMatrixCSC{Tval, Tind}(ldl.m, ldl.n, ldl.colptr_B, ldl.rowval_B, ldl.nzval_B);
    L = [LA;LB];
    L = hcat(L, SparseArrays.sparse(Tind[], Tind[], Tval[], ldl.n + ldl.m, ldl.m));
    
    #get the full permutation vector
    P = [ldl.col; collect(ldl.n+1:(ldl.n+ldl.m))];
    #turn subset from an array into a true subset
    S = Set(subset);
    #get the permutation vector inside subsets
    idxSet = Tind[];
    for ii in 1:1:(ldl.n+ldl.m)
        if(in(P[ii], S))
            push!(idxSet, ii);
        end
    end
    #find the permutation indices in the subset
    PSub = Tind[];
    PInSet = ldl.col[idxSet[1:numInternalNode]];
    for ii in 1:1:numInternalNode
        push!(PSub, subsetMap[PInSet[ii]]);
    end

    #taking only the subset
    LSub = L[idxSet, idxSet];
    
    LSub = LSub[:, 1:numInternalNode];
    LASub = LSub[1:numInternalNode, :];
    LBSub = LSub[numInternalNode + 1:end, :];
    ldlSub = Laplacians.LDL(numInternalNode, numPort, 
                            numInternalNode + 1,
                            PSub, 
                            #LASub section
                            size(LASub.rowval, 1) + 1, 
                            LASub.colptr, 
                            LASub.rowval,
                            LASub.nzval,
                            #LBSub section
                            size(LBSub.rowval, 1) + 1, 
                            LBSub.colptr, 
                            LBSub.rowval,
                            LBSub.nzval, 
                            #diagonal section
                            ldl.diagA[idxSet[1:numInternalNode]]);
    #creating the corresponding subset Schur Complement
    #here using subset or idxSet is the same since we do not permute ports
    schurCSub = schurC[subset[numInternalNode + 1: end] .- ldl.n, subset[numInternalNode + 1: end] .- ldl.n];
    return ldlSub, schurCSub
end

#===========================================

  Alternate solver approach

===========================================#


"""
    L = ldli2Chol(ldli)
This produces a matrix L so that L L^T approximate the original Laplacians.
It is not quite a Cholesky factor, because it is off by a perm
(and the all-1s vector orthogonality.
"""
function ldli2Chol(ldli)
    n = length(ldli.colptr)
    m = n + length(ldli.fval)
    li = zeros(Int,m)
    lj = zeros(Int,m)
    lv = zeros(Float64,m)
    lptr = 0

    dhi = zeros(n)
    for i in 1:n
        if ldli.d[i] == 0
            dhi[i] = 1.0
        else
            dhi[i] = sqrt(ldli.d[i])
        end
    end

    scales = ones(n)
    for ii in 1:(n-1)
        i = ldli.col[ii]
        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-1
        scales[i] = prod(1.0 .- ldli.fval[j0:(j1-1)])
    end

    for ii in 1:(n-1)
        i = ldli.col[ii]
        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-1
        scale = scales[i] / dhi[i]

        scj = 1
        for jj in j0:(j1-1)
            j = ldli.rowval[jj]
            f = ldli.fval[jj]

            lptr += 1
            li[lptr] = i
            lj[lptr] = j
            lv[lptr] = -f*scj/scale


            scj = scj*(1-f)
        end
        j = ldli.rowval[j1]

        lptr += 1
        li[lptr] = i
        lj[lptr] = j
        lv[lptr] = -dhi[i]

        lptr += 1
        li[lptr] = i
        lj[lptr] = i
        lv[lptr] = 1/scale

    end

    for i in 1:n
        if ldli.d[i] == 0
            lptr += 1
            li[lptr] = i
            lj[lptr] = i
            lv[lptr] = 1.0
        end
    end

    return sparse(li,lj,lv,n,n)
    #return li, lj, lv
end

function LDLsolver(L::SparseMatrixCSC, b::Array)
    y = x6 = L \ (L' \ b)
    return y .- mean(y)
end


"""
This variation of approxChol creates a cholesky factor to do the elimination.
It has not yet been optimized, and does not yet make the cholesky factor lower triangular
"""
function approxchol_lapChol(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[]) where {Tv,Ti}

    tol_ =tol
    maxits_ =maxits
    maxtime_ =maxtime
    verbose_ =verbose
    pcgIts_ =pcgIts

    t1 = time()
    llmat = LLmatp(a)

    ldli = approxChol(llmat)

    chL = ldli2Chol(ldli)

    if verbose
      println("Factorization time: ", time()-t1)
      println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
    end

    F(b) = LDLsolver(chL, b)

    la = lap(a)

    if verbose
        println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
    end


    f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b .- mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)

    return f
end




#=============================================================

ApproxCholPQ
It only implements pop, increment key, and decrement key.
All nodes with degrees 1 through n appear in their own doubly-linked lists.
Nodes of higher degrees are bundled together.

=============================================================#

function keyMap(x, n)
    return x <= n ? x : n + div(x,n)
end

function ApproxCholPQ(a::Vector{Tind}) where Tind

    n = length(a)
    elems = Array{ApproxCholPQElem{Tind}}(undef, n)
    lists = zeros(Tind, 2*n+1)
    minlist = one(n)

    for i in 1:length(a)
        key = a[i]
        head = lists[key]

        if head > zero(Tind)
            elems[i] = ApproxCholPQElem{Tind}(zero(Tind), head, key)

            elems[head] = ApproxCholPQElem{Tind}(i, elems[head].next, elems[head].key)
        else
            elems[i] = ApproxCholPQElem{Tind}(zero(Tind), zero(Tind), key)

        end

        lists[key] = i
    end

    return ApproxCholPQ(elems, lists, minlist, n, n)
end
#constructor given size
function ApproxCholPQ(a::Vector{Tind}, num::Tind) where Tind

    n = num;
    elems = Array{ApproxCholPQElem{Tind}}(undef, n)
    lists = zeros(Tind, 2*n+1)
    minlist = one(n)

    @assert n<=length(a)
    for i in 1:n
        key = a[i]
        head = lists[key]

        if head > zero(Tind)
            elems[i] = ApproxCholPQElem{Tind}(zero(Tind), head, key)

            elems[head] = ApproxCholPQElem{Tind}(i, elems[head].next, elems[head].key)
        else
            elems[i] = ApproxCholPQElem{Tind}(zero(Tind), zero(Tind), key)

        end

        lists[key] = i
    end

    return ApproxCholPQ(elems, lists, minlist, n, n)
end

function approxCholPQPop!(pq::ApproxCholPQ{Tind}) where Tind
    if pq.nitems == 0
        error("ApproxPQ is empty")
    end
    while pq.lists[pq.minlist] == 0
        pq.minlist = pq.minlist + 1
    end
    i = pq.lists[pq.minlist]
    next = pq.elems[i].next


    pq.lists[pq.minlist] = next
    if next > 0
        pq.elems[next] = ApproxCholPQElem(zero(Tind), pq.elems[next].next, pq.elems[next].key)
    end

    pq.nitems -= 1

    return i
end

#the function to print an entire priority queue
function PrintApproxCholPQ(pq::ApproxCholPQ{Tind}) where Tind
    println("minlist: $(pq.minlist), nitems: $(pq.nitems), n: $(pq.n)");
    for i in 1:2*pq.n+1
        print("list of deg $(i): ")
        ptr = pq.lists[i]
        while ptr>zero(Tind)
            print("$(ptr)->")
            ptr = pq.elems[ptr].next;
        end
        println("")
    end
end


function approxCholPQMove!(pq::ApproxCholPQ{Tind}, i, newkey, oldlist, newlist) where Tind

    prev = pq.elems[i].prev
    next = pq.elems[i].next

    # remove i from its old list
    if next > zero(Tind)
        pq.elems[next] = ApproxCholPQElem{Tind}(prev, pq.elems[next].next, pq.elems[next].key)
    end
    if prev > zero(Tind)
        pq.elems[prev] = ApproxCholPQElem{Tind}(pq.elems[prev].prev, next, pq.elems[prev].key)

    else
        pq.lists[oldlist] = next
    end

    # insert i into its new list
    head = pq.lists[newlist]
    if head > 0
        pq.elems[head] = ApproxCholPQElem{Tind}(i, pq.elems[head].next, pq.elems[head].key)
    end
    pq.lists[newlist] = i

    pq.elems[i] = ApproxCholPQElem{Tind}(zero(Tind), head, newkey)

    return nothing
end

"""
    Decrement the key of element i
    This could crash if i exceeds the maxkey
"""
function approxCholPQDec!(pq::ApproxCholPQ{Tind}, i) where Tind
    
    if(i>=pq.n)
        return nothing;
    end
    oldlist = keyMap(pq.elems[i].key, pq.n)
    newlist = keyMap(pq.elems[i].key - one(Tind), pq.n)

    if newlist != oldlist

        approxCholPQMove!(pq, i, pq.elems[i].key - one(Tind), oldlist, newlist)

        if newlist < pq.minlist
            pq.minlist = newlist
        end

    else
        pq.elems[i] = ApproxCholPQElem{Tind}(pq.elems[i].prev, pq.elems[i].next, pq.elems[i].key - one(Tind))
    end


    return nothing
end

"""
    Increment the key of element i
    This could crash if i exceeds the maxkey
"""
function approxCholPQInc!(pq::ApproxCholPQ{Tind}, i) where Tind

    if(i>=pq.n)
        return nothing;
    end

    oldlist = keyMap(pq.elems[i].key, pq.n)
    newlist = keyMap(pq.elems[i].key + one(Tind), pq.n)

    if newlist != oldlist

        approxCholPQMove!(pq, i, pq.elems[i].key + one(Tind), oldlist, newlist)

    else
        pq.elems[i] = ApproxCholPQElem{Tind}(pq.elems[i].prev, pq.elems[i].next, pq.elems[i].key + one(Tind))
    end

    return nothing
end
