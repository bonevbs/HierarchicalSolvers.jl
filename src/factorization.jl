### routines to generate the factorization

using Infiltrator

## Symbolic factorization routine
# returns a reordered nessted dissection tree as well as a nested dissection tree containing the local indices
symfact!(nd::NestedDissection) = _symfact!(nd)
function _symfact!(nd::NestedDissection)
  if isleaf(nd) 
    nd_loc = NDNode(Vector{Int}(), Vector{Int}())
  else
    if !isnothing(nd.left)
      nd.left, left_loc = _symfact!(nd.left)
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
      nd.right, right_loc = _symfact!(nd.right)
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
  return nd, nd_loc
end

## Factorization routine
function factor(A::SparseMatrixCSC{T}, nd::NestedDissection) where T
  nd, nd_loc = symfact!(nd)
  nd_loc.int = collect(1:length(nd.bnd))
  nd_loc.bnd = Vector{Int}()
  F = _factor_branch(A, nd, nd_loc)
end

# recursive definition of the internal factorization routine
function _factor_branch(A::SparseMatrixCSC{T}, nd::NestedDissection, nd_loc::NestedDissection) where T
  if isleaf(nd)
    F = _factor_leaf(A, nd.int, nd.bnd, nd_loc.int, nd_loc.bnd)
  elseif isbranch(nd)
    Fl = _factor_branch(A, nd.left, nd_loc.left)
    Fr = _factor_branch(A, nd.right, nd_loc.right)

    int1 = nd.left.bnd[nd_loc.left.int]; bnd1 = nd.left.bnd[nd_loc.left.bnd];
    int2 = nd.right.bnd[nd_loc.right.int]; bnd2 = nd.right.bnd[nd_loc.right.bnd]; 
    ni1 = length(nd_loc.left.int); nb1 = length(nd_loc.left.bnd)
    ni2 = length(nd_loc.right.int); nb2 = length(nd_loc.right.bnd)

    Aii = [Fl.S[1:ni1, 1:ni1] Matrix(A[int1, int2]); Matrix(A[int2, int1]) Fr.S[1:ni2, 1:ni2]];
    Aib = [Fl.S[1:ni1, ni1+1:end] Matrix(A[int1, bnd2]); Matrix(A[int2, bnd1]) Fr.S[1:ni2, ni2+1:end]];
    Abi = [Fl.S[ni1+1:end, 1:ni1] Matrix(A[bnd1, int2]); Matrix(A[bnd2, int1]) Fr.S[ni2+1:end, 1:ni2]];
    Abb = [Fl.S[ni1+1:end, ni1+1:end] Matrix(A[bnd1, bnd2]); Matrix(A[bnd2, bnd1]) Fr.S[ni2+1:end, ni2+1:end]];

    L = Abi / Aii
    R = Aii \ Aib
    S = Abb - Abi * R
    perm = [nd_loc.int; nd_loc.bnd]

    F = FactorNode(Aii, S[perm,perm], L, R, nd.int, nd.bnd, nd_loc.int, nd_loc.bnd, Fl, Fr) # remove local branch storage
  else
    throw(ErrorException("Expected nested dissection to be a binary tree. Found a node with only one child."))  
  end
  return F
end

function _factor_leaf(A::SparseMatrixCSC{T, Int}, int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int}) where T
  D = A[int, int]
  L = Matrix(A[bnd, int]) / D # converts right/left-hand side to dense first
  R = D \ Matrix(A[int, bnd])
  S = A[bnd, bnd] - A[bnd, int] * R
  #S = A[bnd[[int_loc; bnd_loc]],bnd[int_loc; bnd_loc]] - A[bnd[int_loc; bnd_loc], int] * R[:, [int_loc; bnd_loc]]
  perm = [int_loc; bnd_loc]
  F = FactorNode(D, S[perm,perm], L, R, int, bnd, int_loc, bnd_loc)
  return F
end
