### routines to generate the factorization
# Written by Boris Bonev, Feb. 2021

using Infiltrator

## Factorization routine
function factor(A::SparseMatrixCSC{T}, nd::NestedDissection, opts::SolverOptions=SolverOptions(T);  args...) where T
  opts = copy(opts; args...)
  chkopts!(opts)
  nd_loc = symfact!(nd)
  nd_loc.int = collect(1:length(nd.bnd))
  nd_loc.bnd = Vector{Int}()
  F = _factor(A, nd, nd_loc, 1; opts.swlevel, opts.atol, opts.rtol)
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
function _factor(A::AbstractMatrix{T}, nd::NestedDissection, nd_loc::NestedDissection, level::Int; swlevel::Int, atol::Float64, rtol::Float64) where T
  if isleaf(nd)
    F = _factor_leaf(A, nd.int, nd.bnd, nd_loc.int, nd_loc.bnd; compress=level≤swlevel, atol=atol, rtol=rtol)
  elseif isbranch(nd)
    Fl = _factor(A, nd.left, nd_loc.left, level+1; swlevel, atol, rtol)
    Fr = _factor(A, nd.right, nd_loc.right, level+1; swlevel, atol, rtol)
    F = _factor_branch(A, Fl, Fr, nd, nd_loc; compress=level≤swlevel, atol=atol, rtol=rtol)
  else
    throw(ErrorException("Expected nested dissection to be a binary tree. Found a node with only one child."))  
  end
  return F
end

function _factor_leaf(A::AbstractMatrix{T}, int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int}; compress::Bool, atol::Float64, rtol::Float64) where T
  D = A[int, int]
  L = Matrix(A[bnd, int]) / D # converts right/left-hand side to dense first
  R = D \ Matrix(A[int, bnd])
  S = A[bnd, bnd] - A[bnd, int] * R
  #S = A[bnd[[int_loc; bnd_loc]],bnd[int_loc; bnd_loc]] - A[bnd[int_loc; bnd_loc], int] * R[:, [int_loc; bnd_loc]]
  perm = [int_loc; bnd_loc]
  F = FactorNode(D, S[perm,perm], L, R, int, bnd, int_loc, bnd_loc)
  return F
end

# this is where the magic happens
function _factor_branch(A::AbstractMatrix{T}, F1::FactorNode{T}, F2::FactorNode{T}, nd::NestedDissection, nd_loc::NestedDissection; compress::Bool, atol::Float64, rtol::Float64) where T
  int1 = nd.left.bnd[nd_loc.left.int]; bnd1 = nd.left.bnd[nd_loc.left.bnd];
  int2 = nd.right.bnd[nd_loc.right.int]; bnd2 = nd.right.bnd[nd_loc.right.bnd]; 
  ni1 = length(nd_loc.left.int); nb1 = length(nd_loc.left.bnd)
  ni2 = length(nd_loc.right.int); nb2 = length(nd_loc.right.bnd)

  S1 = F1.S; S2 = F2.S

  # TODO: Split this into two parts: one for Hss, one for normal matrices
  if typeof(F1.S) <: HssMatrix || typeof(F2.S) <: HssMatrix
    # TODO: check that the blocking is actually
    rcl1, ccl1 = cluster(F1.S)
    rcl2, ccl2 = cluster(F2.S)
    # extreact generators of children Schur complements
    Uint1, Vint1 = generators(F1.S.A11); Uint1 .= Uint1*S1.B12
    Ubnd1, Vbnd1 = generators(F1.S.A22); Ubnd1 .= Ubnd1*S1.B21
    Uint2, Vint2 = generators(F2.S.A11); Uint2 .= Uint2*S2.B12
    Ubnd2, Vbnd2 = generators(F2.S.A22); Ubnd2 .= Ubnd2*S2.B12
    # form the blocks
    Aii = BlockMatrix(F1.S.A11, compress(A[int1, int2], rcl1, ccl2), compress(A[int2, int1], rcl2, ccl1), F2.S.A11) # check hssranks of the offdiagonal guys
    Aib = BlockMatrix(LowRankMatrix(Uint1, Vbnd1), A[int1, bnd2], A[int2, bnd1], LowRankMatrix(Uint2, Vbnd2))
    Abi = BlockMatrix(LowRankMatrix(Ubnd1, Vint1), A[bnd1, int2], A[bnd2, int1], LowRankMatrix(Ubnd2, Vint2))
    Abb = BlockMatrix(F1.S.A22, A[bnd1, bnd2], A[bnd2, bnd1], F2.S.A22)
  else # save everything densely
    Aii = BlockMatrix(view(F1.S, 1:ni1, 1:ni1), view(A, int1, int2), view(A, int2, int1), view(F2.S, 1:ni2, 1:ni2))
    Aib = BlockMatrix(view(F1.S, 1:ni1, ni1+1:ni1+nb1), view(A, int1, bnd2), view(A, int2, bnd1), view(F2.S, 1:ni2, ni2+1:ni2+nb2))
    Abi = BlockMatrix(view(F1.S, ni1+1:ni1+nb1, 1:ni1), view(A, bnd1, int2), view(A, bnd2, int1), view(F2.S, ni2+1:ni2+nb2, 1:ni2))
    Abb = BlockMatrix(view(F1.S, ni1+1:ni1+nb1, ni1+1:ni1+nb1), view(A, bnd1, bnd2), view(A, bnd2, bnd1), view(F2.S, ni2+1:ni2+nb2, ni2+1:ni2+nb2))
  end

  # Form the Factorization by forming the Gauss transforms and the Schur complements
  if compress
    # build operators
    Lmul = (y, _, x) ->  y = Abi*(Aii\x)
    Lmulc = (y, _, x) ->  y = ((x'*Abi)/Aii)'
    Lop = LinearOperator{T}(nb1+nb2, ni1+ni2, Lmul, Lmulc, nothing);
    Rmul = (y, _, x) ->  y = Aii\(Aib*x)
    Rmulc = (y, _, x) ->  y = ((x'/Aii)*Aib)'
    Rop = LinearOperator{T}(ni1+ni2, nb1+nb2, Rmul, Rmulc, nothing);
    #@infiltrate
    # perform the sampling
    F = pqrfact(Lop, sketch=:randn, atol=atol, rtol=rtol)
    L = LowRankMatrix(F.Q, collect(F.R[:,F.p]'))
    F = pqrfact(Rop, sketch=:randn, atol=atol, rtol=rtol)
    R = LowRankMatrix(F.Q, collect(F.R[:,F.p]'))

    S = Matrix(Abb) - Matrix(Abi) * Matrix(R)
    perm = [nd_loc.int; nd_loc.bnd];
    F = FactorNode(Matrix(Aii), S[perm,perm], L, R, nd.int, nd.bnd, nd_loc.int, nd_loc.bnd, F1, F2) # remove local branch storage
  else
    L = Abi/Aii
    R = Aii\Aib
    #println(size(R))
    #println(size(Abb))
    S = Abb - Abi*R
    perm = [nd_loc.int; nd_loc.bnd];
    F = FactorNode(Matrix(Aii), S[perm,perm], Matrix(L), Matrix(R), nd.int, nd.bnd, nd_loc.int, nd_loc.bnd, F1, F2) # remove local branch storage
  end

  return F

  # to compress or not to compress
end