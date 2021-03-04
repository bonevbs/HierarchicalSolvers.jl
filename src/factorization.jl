### routines to generate the factorization
# Written by Boris Bonev, Feb. 2021

using Infiltrator

## Factorization routine
function factor(A::SparseMatrixCSC{T}, nd::NestedDissection, opts::SolverOptions=SolverOptions(T);  args...) where T
  opts = copy(opts; args...)
  chkopts!(opts)
  opts.swlevel < 0 ? swlevel = max(depth(nd) + opts.swlevel, 0) : swlevel = opts.swlevel
  nd_loc = symfact!(nd)
  nd_loc.int = collect(1:length(nd.bnd))
  nd_loc.bnd = Vector{Int}()
  F = _factor(A, nd, nd_loc, 1; swlevel, opts.atol, opts.rtol, opts.leafsize)
  return F
end

## Symbolic factorization routine
# returns a reordered nessted dissection tree as well as a nested dissection tree containing the local indices
symfact!(nd::NestedDissection) = _symfact!(nd, 1)
function _symfact!(nd::NestedDissection, level)
  if isleaf(nd)
    nd_loc = NDNode(Vector{Int}(), Vector{Int}())
  else
    if !isnothing(nd.left)
      left_loc = _symfact!(nd.left, level+1)
      # figure out where the children indices go in the parent
      left_loc.int = findall(in(nd.int), nd.left.bnd)
      left_loc.bnd = findall(in(nd.bnd), nd.left.bnd)
      # reordered indices for the current node
      intl = nd.left.bnd[left_loc.int];
      bndl = nd.left.bnd[left_loc.bnd];
    else
      intl = Vector{Int}()
      bndl = Vector{Int}()
      left_loc = nothing
    end
    if !isnothing(nd.right)
      right_loc = _symfact!(nd.right, level+1)
      right_loc.int = findall(in(nd.int), nd.right.bnd)
      right_loc.bnd = findall(in(nd.bnd), nd.right.bnd)
      intr = nd.right.bnd[right_loc.int];
      bndr = nd.right.bnd[right_loc.bnd];
    else
      intr = Vector{Int}()
      bndr = Vector{Int}()
      right_loc = nothing
    end
    # check that we really got all the degrees of freedom
    nd.int = [intl; intr]
    nd.bnd = [bndl; bndr]
    nd_loc = NDNode(Vector{Int}(), Vector{Int}(), left_loc, right_loc)
  end
  return nd_loc
end

# recursive definition of the internal factorization routine
function _factor(A::AbstractMatrix{T}, nd::NestedDissection, nd_loc::NestedDissection, level::Int; swlevel::Int, atol::Float64, rtol::Float64, leafsize::Int) where T
  if isleaf(nd)
    F = _factor_leaf(A, nd.int, nd.bnd, nd_loc.int, nd_loc.bnd, level≤swlevel, atol, rtol, leafsize)
  elseif isbranch(nd)
    Fl = _factor(A, nd.left, nd_loc.left, level+1; swlevel, atol, rtol, leafsize)
    Fr = _factor(A, nd.right, nd_loc.right, level+1; swlevel, atol, rtol, leafsize)
    F = _factor_branch(A, Fl, Fr, nd, nd_loc, level≤swlevel, atol, rtol, leafsize)
  else
    throw(ErrorException("Expected nested dissection to be a binary tree. Found a node with only one child."))  
  end
  return F
end

function _factor_leaf(A::AbstractMatrix{T}, int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int}, cmpflag::Bool, atol::Float64, rtol::Float64, leafsize::Int) where T
  D = A[int, int]
  L = A[bnd, int] / D
  R = D \ A[int, bnd]
  S = zeros(T, length(bnd), length(bnd) )
  perm = [int_loc; bnd_loc]
  S[invperm(perm), invperm(perm)] .= A[bnd, bnd] .- A[bnd, int] * R
  if cmpflag
    cl = bisection_cluster((length(int_loc), length(bnd)); leafsize)
    hssS = compress(S, cl, cl, atol=atol, rtol=rtol)
    F = FactorNode(D, hssS, L, R, int, bnd, int_loc, bnd_loc)
  else
    F = FactorNode(D, S, L, R, int, bnd, int_loc, bnd_loc)
  end
  return F
end

# this is where the magic happens
function _factor_branch(A::AbstractMatrix{T}, Fl::FactorNode{T}, Fr::FactorNode{T}, nd::NestedDissection, nd_loc::NestedDissection, cmpflag::Bool, atol::Float64, rtol::Float64, leafsize::Int) where T
  int1 = nd.left.bnd[nd_loc.left.int]; bnd1 = nd.left.bnd[nd_loc.left.bnd];
  int2 = nd.right.bnd[nd_loc.right.int]; bnd2 = nd.right.bnd[nd_loc.right.bnd]; 

  Aii, Aib, Abi, Abb = _assemble_blocks(A, Fl.S, Fr.S, int1, int2, bnd1, bnd2; atol, rtol)

  # Form the Factorization by forming the Gauss transforms and the Schur complements
  if cmpflag
    # build operators
    Lmul = (y, _, x) ->  y = Abi*(Aii\x); Lmulc = (y, _, x) ->  y = ((x'*Abi)/Aii)'
    Lop = LinearOperator{T}(size(Abi)..., Lmul, Lmulc, nothing);
    Rmul = (y, _, x) ->  y = Aii\(Aib*x); Rmulc = (y, _, x) ->  y = ((x'/Aii)*Aib)'
    Rop = LinearOperator{T}(size(Aib)..., Rmul, Rmulc, nothing);
    # perform the sampling # TODO: replace this with c_tol
    F = pqrfact(Lop, sketch=:randn, atol=0.5*atol, rtol=0.5*rtol)
    L = LowRankMatrix(F.Q, collect(F.R[:,invperm(F.p)]'))
    F = pqrfact(Rop, sketch=:randn, atol=0.5*atol, rtol=0.5*rtol)
    R = LowRankMatrix(F.Q, collect(F.R[:,invperm(F.p)]'))

    #cl = bisection_cluster((length(int_loc), length(bnd)); leafsize)
    S = Matrix(Abb) - Matrix(Abi*R)
    perm = [nd_loc.int; nd_loc.bnd];
    F = FactorNode(Matrix(Aii), S[perm,perm], L, R, nd.int, nd.bnd, nd_loc.int, nd_loc.bnd, Fl, Fr) # remove local branch storage
  else
    L = Abi/Aii
    R = Aii\Aib
    S = Abb - Abi*R
    perm = [nd_loc.int; nd_loc.bnd];
    F = FactorNode(Matrix(Aii), S[perm,perm], Matrix(L), Matrix(R), nd.int, nd.bnd, nd_loc.int, nd_loc.bnd, Fl, Fr) # remove local branch storage
  end

  return F

  # to compress or not to compress
end

function _assemble_blocks(A::AbstractMatrix{T}, S1::AbstractMatrix{T}, S2::AbstractMatrix{T}, int1::Vector{Int}, int2::Vector{Int}, bnd1::Vector{Int}, bnd2::Vector{Int}; atol::Float64, rtol::Float64) where T
  ni1 = length(int1); nb1 = length(bnd1)
  ni2 = length(int2); nb2 = length(bnd2)
  Aii = BlockMatrix(S1[1:ni1, 1:ni1], A[int1, int2], A[int2, int1], S2[1:ni2, 1:ni2])
  Aib = BlockMatrix(S1[1:ni1, ni1+1:ni1+nb1], A[int1, bnd2], A[int2, bnd1], S2[1:ni2, ni2+1:ni2+nb2])
  Abi = BlockMatrix(S1[ni1+1:ni1+nb1, 1:ni1], A[bnd1, int2], A[bnd2, int1], S2[ni2+1:ni2+nb2, 1:ni2])
  Abb = BlockMatrix(S1[ni1+1:ni1+nb1, ni1+1:ni1+nb1], A[bnd1, bnd2], A[bnd2, bnd1], S2[ni2+1:ni2+nb2, ni2+1:ni2+nb2])
  return Aii, Aib, Abi, Abb
end

function _assemble_blocks(A::AbstractMatrix{T}, S1::HssMatrix{T}, S2::HssMatrix{T}, int1::Vector{Int}, int2::Vector{Int}, bnd1::Vector{Int}, bnd2::Vector{Int}; atol::Float64, rtol::Float64) where T
  ni1 = length(int1); ni2 = length(int2); nb1 = length(bnd1); nb2 = length(bnd2)
  # TODO: Split this into two parts: one for Hss, one for normal matrices
  # TODO: move this into it's own block for performance
  # TODO: check that the blocking is actually
  rcl1, ccl1 = cluster(S1.A11); rcl2, ccl2 = cluster(S2.A11)
  # extract generators of children Schur complements
  Uint1, Vint1 = generators(S1.A11); Uint1 = Uint1*S1.B12
  Uint2, Vint2 = generators(S2.A11); Uint2 = Uint2*S2.B12
  Ubnd1, Vbnd1 = generators(S1.A22); Ubnd1 = Ubnd1*S1.B21
  Ubnd2, Vbnd2 = generators(S2.A22); Ubnd2 = Ubnd2*S2.B21
  # form the blocks
  Aii = BlockMatrix(S1.A11, hss(A[int1, int2], rcl1, ccl2, atol=atol, rtol=rtol), hss(A[int2, int1], rcl2, ccl1, atol=atol, rtol=rtol), S2.A11) # check hssranks of the offdiagonal guys
  Aib = BlockMatrix(LowRankMatrix(Uint1, Vbnd1), A[int1, bnd2], A[int2, bnd1], LowRankMatrix(Uint2, Vbnd2))
  Abi = BlockMatrix(LowRankMatrix(Ubnd1, Vint1), A[bnd1, int2], A[bnd2, int1], LowRankMatrix(Ubnd2, Vint2))
  Abb = BlockMatrix(S1.A22, A[bnd1, bnd2], A[bnd2, bnd1], S2.A22)
  return Aii, Aib, Abi, Abb
end