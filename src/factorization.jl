### routines to generate the factorization
# Written by Boris Bonev, Feb. 2021

## Factorization routine
function factor(A::SparseMatrixCSC{T}, nd::NestedDissection, nd_loc::NestedDissection, opts::SolverOptions=SolverOptions();  args...) where T
  opts = copy(opts; args...)
  chkopts!(opts)
  opts.swlevel < 0 ? swlevel = max(depth(nd) + opts.swlevel, 0) : swlevel = opts.swlevel
  F = _factor(A, nd, nd_loc, 1; swlevel, opts.swsize, opts.atol, opts.rtol, opts.leafsize, opts.kest, opts.stepsize, opts.verbose)
  return F
end

# recursive definition of the internal factorization routine
function _factor(A::AbstractMatrix{T}, nd::NestedDissection, nd_loc::NestedDissection, level::Int; swlevel::Int, swsize::Int, atol::Float64, rtol::Float64, leafsize::Int, kest::Int, stepsize::Int, verbose::Bool) where T
  compression_flag = (level≤swlevel) && (length(nd.bnd)≥swsize) # determine whether compression is needed for this level
  if isleaf(nd)
    verbose && println("Factoring leafnode at level ", level)
    return _factor_leaf(A, nd, nd_loc, Val(compression_flag); atol, rtol, leafsize)
  elseif isbranch(nd)
    Fl = _factor(A, nd.left, nd_loc.left, level+1; swlevel, swsize, atol, rtol, leafsize, kest, stepsize, verbose)
    Fr = _factor(A, nd.right, nd_loc.right, level+1; swlevel, swsize, atol, rtol, leafsize, kest, stepsize, verbose)
    verbose && println("Factoring branchnode at level ", level)
    return _factor_branch(A, Fl, Fr, nd, nd_loc, Val(compression_flag); atol, rtol, leafsize, kest, stepsize, verbose)
  else
    throw(ErrorException("Expected nested dissection to be a binary tree. Found a node with only one child."))  
  end
end

# factor leaf node without compressing it
function _factor_leaf(A::AbstractMatrix{T}, nd::NestedDissection, nd_loc::NestedDissection, ::Val{false}; args...) where T
  @timeit to "factor diagonal block (leaves)" begin
    int = nd.int; bnd = nd.bnd;
    int_loc = nd_loc.int; bnd_loc = nd_loc.bnd
    D = Matrix(view(A, int, int))
  end
  @timeit to "Gauss transforms (leaves)" begin
    Abi = Matrix(view(A, bnd, int))
    L = Abi / D
    R = D \ Matrix(view(A,int, bnd))
  end
  @timeit to "Schur complement (leaves)" begin
    perm = [int_loc; bnd_loc]
    S = Matrix(view(A,bnd, bnd)) .- Abi * R
  end
  return FactorNode(D, S[perm, perm], L, R, int, bnd, int_loc, bnd_loc)
end

# factor leaf node and compress to HSS form (almost never gets called unless compression starts at the bottom level)
function _factor_leaf(A::AbstractMatrix{T}, nd::NestedDissection, nd_loc::NestedDissection, ::Val{true}; atol::Float64, rtol::Float64, leafsize::Int) where T
  @timeit to "factor diagonal block (leaves)" begin  
    int = nd.int; bnd = nd.bnd;
    int_loc = nd_loc.int; bnd_loc = nd_loc.bnd
    D = Matrix(view(A, int, int))
  end
  @timeit to "Gauss transforms (leaves)" begin
    Abi = Matrix(view(A, bnd, int))
    L = Matrix(view(A, bnd, int)) / D
    R = D \ Matrix(view(A,int, bnd))
  end
  @timeit to "Schur complement (leaves)" begin
    perm = [int_loc; bnd_loc]
    S = Matrix(view(A,bnd, bnd)) .- Abi * R
    cl = bisection_cluster((length(int_loc), length(bnd)); leafsize)
    hssS = compress(S[perm, perm], cl, cl; atol=atol, rtol=rtol)
  end
  return FactorNode(D, hssS, L, R, int, bnd, int_loc, bnd_loc)
end

# factor branch node without compressing it
function _factor_branch(A::AbstractMatrix{T}, Fl::FactorNode{T}, Fr::FactorNode{T}, nd::NestedDissection, nd_loc::NestedDissection, ::Val{false}; atol::Float64, rtol::Float64, leafsize::Int, kest::Int, stepsize::Int, verbose::Bool) where T
  int1 = nd.left.bnd[nd_loc.left.int]; bnd1 = nd.left.bnd[nd_loc.left.bnd];
  int2 = nd.right.bnd[nd_loc.right.int]; bnd2 = nd.right.bnd[nd_loc.right.bnd]; 
  int_loc = nd_loc.int; bnd_loc = nd_loc.bnd
  
  @timeit to "assemble blocks" begin
    Aii, Aib, Abi, Abb = _assemble_blocks(A, Fl.S, Fr.S, int1, int2, bnd1, bnd2; atol, rtol)
  end

  @timeit to "factor diagonal block" begin
    Aii = blockfactor(Aii; atol, rtol)
  end
  @timeit to "Gauss transforms" begin
    L = blockrdiv(Abi, Aii)
    R = blockldiv(Aii, Aib)
  end
  @timeit to "Schur complement" begin
    S = Abb - Abi*R
    perm = [nd_loc.int; nd_loc.bnd];
  end
  return FactorNode(Aii, S[perm,perm], L, R, nd.int, nd.bnd, nd_loc.int, nd_loc.bnd, Fl, Fr)
end

# factor node matrix free and compress it
function _factor_branch(A::AbstractMatrix{T}, Fl::FactorNode{T}, Fr::FactorNode{T}, nd::NestedDissection, nd_loc::NestedDissection, ::Val{true}; atol::Float64, rtol::Float64, leafsize::Int, kest::Int, stepsize::Int, verbose::Bool) where T
  int1 = contigious(nd.left.bnd[nd_loc.left.int]); bnd1 = contigious(nd.left.bnd[nd_loc.left.bnd]);
  int2 = contigious(nd.right.bnd[nd_loc.right.int]); bnd2 = contigious(nd.right.bnd[nd_loc.right.bnd]); 
  # int1 = nd.left.bnd[nd_loc.left.int]; bnd1 = nd.left.bnd[nd_loc.left.bnd];
  # int2 = nd.right.bnd[nd_loc.right.int]; bnd2 = nd.right.bnd[nd_loc.right.bnd]; 
  int_loc = nd_loc.int; bnd_loc = nd_loc.bnd

  # make sure clusters are matching
  if ishss(Fl.S) && ishss(Fr.S)
    S1, S2 = _equilibrate_clusters(Fl.S, Fr.S)
  else
    S1, S2 = Fl.S, Fr.S
  end
  @timeit to "assemble blocks (compressed)" begin
    Aii, Aib, Abi, Abb = _assemble_blocks(A, S1, S2, int1, int2, bnd1, bnd2; atol, rtol)
  end

  # block-factorization
  @timeit to "factor diagonal block (compressed)" begin
    D = blockfactor(Aii; atol, rtol)
  end

  # use randomized compression to get the low-rank representation
  # TODO: replace this with c_tol
  @timeit to "Gauss transforms (compressed)" begin
    # build operators
    Lop, Rop = _gauss_transforms(D, Aib, Abi)

    # this is the final part which is type unstable!!!!!
    qrfL = pqrfact(Lop, sketch=:randn, atol=0.5*atol, rtol=0.5*rtol)
    L = LowRankMatrix(qrfL.Q, collect(qrfL.R[:,invperm(qrfL.p)]'))
    qrfR = pqrfact(Rop, sketch=:randn, atol=0.5*atol, rtol=0.5*rtol)
    R = LowRankMatrix(qrfR.Q, collect(qrfR.R[:,invperm(qrfR.p)]'))
  end

  # use randomized compression to compute the HSS form of the Schur complement
  @timeit to "Schur complement (compressed)" begin
    perm = [nd_loc.int; nd_loc.bnd]
    Smap = _schur_complement(Abb, Abi, R, perm)
    cl = bisection_cluster((length(int_loc), length(int_loc)+length(bnd_loc)); leafsize)
    hssS = randcompress_adaptive(Smap, cl, cl; kest=kest, atol=atol, rtol=rtol, verbose=verbose)
  end
  return FactorNode(D, hssS, L, R, nd.int, nd.bnd, nd_loc.int, nd_loc.bnd, Fl, Fr)
end

# general routine for assembling the new block matrices
function _assemble_blocks(A::AbstractMatrix{T}, S1::AbstractMatrix{T}, S2::AbstractMatrix{T}, int1::AbstractVector{Int}, int2::AbstractVector{Int}, bnd1::AbstractVector{Int}, bnd2::AbstractVector{Int}; args...) where T
  ni1 = length(int1); nb1 = length(bnd1)
  ni2 = length(int2); nb2 = length(bnd2)
  Aii = BlockMatrix(S1[1:ni1, 1:ni1], Matrix(view(A, int1, int2)), Matrix(view(A, int2, int1)), S2[1:ni2, 1:ni2])
  Aib = BlockMatrix(S1[1:ni1, ni1+1:ni1+nb1], Matrix(view(A, int1, bnd2)), Matrix(view(A, int2, bnd1)), S2[1:ni2, ni2+1:ni2+nb2])
  Abi = BlockMatrix(S1[ni1+1:ni1+nb1, 1:ni1], Matrix(view(A, bnd1, int2)), Matrix(view(A, bnd2, int1)), S2[ni2+1:ni2+nb2, 1:ni2])
  Abb = BlockMatrix(S1[ni1+1:ni1+nb1, ni1+1:ni1+nb1], Matrix(view(A, bnd1, bnd2)), Matrix(view(A, bnd2, bnd1)), S2[ni2+1:ni2+nb2, ni2+1:ni2+nb2])
  return Aii, Aib, Abi, Abb
end

# For HSS matrices we want to specialize the routine in order to exploit the pre-determined blocking which exposes interior DOFs
function _assemble_blocks(A::AbstractMatrix{T}, S1::HssMatrix{T}, S2::HssMatrix{T}, int1::AbstractVector{Int}, int2::AbstractVector{Int}, bnd1::AbstractVector{Int}, bnd2::AbstractVector{Int}; atol::Float64, rtol::Float64, verbose=false) where T
  rcl1, ccl1 = cluster(S1.A11); rcl2, ccl2 = cluster(S2.A11)
  # extract generators of children Schur complements
  Uint1, Vint1 = generators(S1.A11); Uint1 = Uint1*S1.B12
  Uint2, Vint2 = generators(S2.A11); Uint2 = Uint2*S2.B12
  Ubnd1, Vbnd1 = generators(S1.A22); Ubnd1 = Ubnd1*S1.B21
  Ubnd2, Vbnd2 = generators(S2.A22); Ubnd2 = Ubnd2*S2.B21

  # form the blocks
  Aii = BlockMatrix(S1.A11, hss(A[int1, int2], rcl1, ccl2; atol=atol, rtol=rtol), hss(A[int2, int1], rcl2, ccl1; atol=atol, rtol=rtol), S2.A11) # check hssranks of the offdiagonal guys
  Aib = BlockMatrix(LowRankMatrix(Uint1, Vbnd1), A[int1, bnd2], A[int2, bnd1], LowRankMatrix(Uint2, Vbnd2))
  Abi = BlockMatrix(LowRankMatrix(Ubnd1, Vint1), A[bnd1, int2], A[bnd2, int1], LowRankMatrix(Ubnd2, Vint2))
  Abb = BlockMatrix(S1.A22, A[bnd1, bnd2], A[bnd2, bnd1], S2.A22)
  return Aii, Aib, Abi, Abb
end

# make sure HSS structures among Schur complements is matching
function _equilibrate_clusters(S1::HssMatrix, S2::HssMatrix)
  rcl1, ccl1 = cluster(S1.A11); rcl2, ccl2 = cluster(S2.A11)
  if !compatible(rcl1, rcl2)
    while !compatible(rcl1, rcl2)
      if depth(rcl1) > depth(rcl2)
        verbose && println("Pruning clusters of node 1")
        rcl1 = prune_leaves!(rcl1); ccl1 = prune_leaves!(ccl1)
        S1.A11 = prune_leaves!(S1.A11)
      elseif depth(rcl1) < depth(rcl2)
        verbose && println("Pruning clusters of node 2")
        rcl2 = prune_leaves!(rcl2); ccl2 = prune_leaves!(ccl2)
        S2.A11 = prune_leaves!(S2.A11)
      else
        verbose && println("Pruning both clusters")
        rcl1 = prune_leaves!(rcl1); ccl1 = prune_leaves!(ccl1)
        rcl2 = prune_leaves!(rcl2); ccl2 = prune_leaves!(ccl2)
        S1.A11 = prune_leaves!(S1.A11)
        S2.A11 = prune_leaves!(S2.A11)
      end
    end
    if isleaf(S1) || isleaf(S2)
      error("One of the Schur complements turned into a leaf. Aborting.")
    end
  end
  return S1, S2
end

function _gauss_transforms(D::BlockFactorization{T}, Aib::BlockMatrix{T}, Abi::BlockMatrix{T}) where T
  Lmul = (y, _, x) ->  y .= Abi*blockldiv!(D, x)
  Lmulc = (y, _, x) ->  y .= blockrdiv!(x'*Abi, D)'
  Rmul = (y, _, x) ->  y .= blockldiv!(D, Aib*x)
  Rmulc = (y, _, x) ->  y .= (blockrdiv!(x', D)*Aib)'
  Lop = LinearOperator{T}(size(Abi)..., Lmul, Lmulc, nothing)
  Rop = LinearOperator{T}(size(Aib)..., Rmul, Rmulc, nothing)
  return Lop, Rop
end

function _schur_complement(Abb::BlockMatrix{T}, Abi::BlockMatrix{T}, R::LowRankMatrix{T}, perm::Vector{Int}) where T
  iperm = invperm(perm)
  U = Abi*R
  Smul = (y, _, x) -> _sample_schur!(y, Abb, U, x, iperm)
  Smulc = (y, _, x) -> _sample_schur!(y, Abb', U', x, iperm)
  Sidx = (i,j) -> _getindex_schur(Abb, U, perm, i, j)
  return LinearMap{T}(size(Abb)..., Smul, Smulc, Sidx)
end

# function for correctly applying the Schur complement
# TODO: this can proabbly be accelerated even further by paying attention to allocation and using mul!
function _sample_schur!(y::AbstractMatrix{T}, A::BlockMatrix{T}, U::LowRankMatrix{T}, x::AbstractMatrix{T}, iperm::Vector{Int}) where T
  #mul!(@view(y[iperm,:]), A, @view(x[iperm,:]), 1., 0.)
  #mul!(@view(y[iperm,:]), B, @view(x[iperm,:]), -1., 1.)
  y[iperm,:] .= A*x[iperm,:] .- U.U*(U.V'*x[iperm,:])
  return y
end

function _getindex_schur(A::AbstractMatrix{T}, U::LowRankMatrix{T}, perm::Vector{Int}, i, j) where T
  ii = perm[i]; jj = perm[j];
  return A[ii, jj] - U.U[ii, :]*U.V[jj,:]'
end
