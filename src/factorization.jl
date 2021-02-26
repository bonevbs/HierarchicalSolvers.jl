### routines to generate the factorization
# Written by Boris Bonev, Feb. 2021

## Factorization routine
function factor(A::SparseMatrixCSC{T}, nd::NestedDissection, opts::SolverOptions=SolverOptions(T);  args...) where T
  opts = copy(opts; args...)
  chkopts!(opts)
  nd_loc = symfact!(nd; swlevel=opts.swlevel)
  nd_loc.int = collect(1:length(nd.bnd))
  nd_loc.bnd = Vector{Int}()
  F = _factor(A, nd, nd_loc, 1; swlevel)
  return F
end

## Symbolic factorization routine
# returns a reordered nessted dissection tree as well as a nested dissection tree containing the local indices
symfact(nd::NestedDissection) = _symfact!(nd, 1)
function _symfact(nd::NestedDissection, level)
  if isleaf(nd)
    nd_loc = NDNode(Vector{Int}(), Vector{Int}())
  else
    if !isnothing(nd.left)
      left_loc = _symfact(nd.left, level+1)
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
function _factor(A::SparseMatrixCSC{T}, nd::NestedDissection, nd_loc::NestedDissection, level::Int; swlevel::Int) where T
  if isleaf(nd)
    F = _factor_leaf(A, nd.int, nd.bnd, nd_loc.int, nd_loc.bnd; level ≤ swlevel)
  elseif isbranch(nd)
    Fl = _factor(A, nd.left, nd_loc.left, level+1; swlevel)
    Fr = _factor(A, nd.right, nd_loc.right, level+1; swlevel)
    F = _factor_branch(A, Fl, Fr, nd, nd_loc; level ≤ swlevel)
  else
    throw(ErrorException("Expected nested dissection to be a binary tree. Found a node with only one child."))  
  end
  return F
end

function _factor_leaf(A::SparseMatrixCSC{T, Int}, int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int}; compress) where T
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
function _factor_branch(A::SparseMatrixCSC{T}, Fl::FactorNode{T}, Fr::FactorNode{T}, nd::NestedDissection, nd_loc::NestedDissection; compress) where T
  int1 = nd.left.bnd[nd_loc.left.int]; bnd1 = nd.left.bnd[nd_loc.left.bnd];
  int2 = nd.right.bnd[nd_loc.right.int]; bnd2 = nd.right.bnd[nd_loc.right.bnd]; 
  ni1 = length(nd_loc.left.int); nb1 = length(nd_loc.left.bnd)
  ni2 = length(nd_loc.right.int); nb2 = length(nd_loc.right.bnd)

  if typeof(Fl.S) <: HssMatrix
    view()
  end

  # TODO: Split this into two parts: one for Hss, one for normal matrices
  if typeof(Fl.S) <: HssMatrix || typeof(Fr.S) <: HssMatrix

  else # save everything densely
    Aii = BlockMatrix(view(Fl.S, 1:ni1, 1:ni1), view(A, int1, int2), view(A, int2, int1), view(Fr.S, 1:ni2, 1:ni2))
    Aib = BlockMatrix(view(Fl.S, 1:ni1, ni1+1:ni1+nb1), view(A, int1, bnd2), view(A, int2, bnd1), view(Fr.S, 1:ni2, ni2+1:ni2+nb2))
    Abi = BlockMatrix(view(Fl.S, ni1+1:ni1+nb1, 1:ni1), view(A, bnd1, int2), view(A, bnd2, int1), view(Fr.S, ni2+1:ni2+nb2, 1:ni2))
    Abb = BlockMatrix(view(Fl.S, ni1+1:ni1+nb1, ni1+1:ni1+nb1), view(A, bnd1, bnd2), view(A, bnd2, bnd1), view(Fr.S, ni2+1:ni2+nb2, ni2+1:ni2+nb2))

    L = Matrix(Abi) / Matrix(Aii)
    R = Matrix(Aii) \ Matrix(Aib)
    S = Matrix(Abb) - Matrix(Abi) * Matrix(R)
    perm = [nd_loc.int; nd_loc.bnd];
    F = FactorNode(Matrix(Aii), S[perm,perm], L, R, nd.int, nd.bnd, nd_loc.int, nd_loc.bnd, Fl, Fr) # remove local branch storage
  end

  # to compress or not to compress
end