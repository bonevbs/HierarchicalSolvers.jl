### routines to generate the factorization

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
function factor(A::SparseMatrixCSC{T}, nd::NestedDissection)
  nd, nd_loc = symfact!(nd)

end

function _factor_branch()
end

function _factor_leaf(A::SparseMatrixCSC{T, Int}, int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int}) where T
  D = A[int, int]
  L = Matrix(A[bnd, int]) / D # converts right/left-hand side to dense first
  R = D \ Matrix(A[int, bnd])
  S = A[bnd,bnd] - A[bnd, int] * R
  # probably unnecessary to save the local indices
  SolverNode(D,S,L,R,int,bnd)
end


# # wrapper to access therecursive definition of the factorization
# function factor(A::AbstractMatrix{T}, nd::NestedDissection)
#   F = _factor_branch(A, nd)
# end

# # recursive definition of the factorization routine
# function _factor_branch(A::AbstractMatrix{T}, nd::NestedDissection, fint::Vector{Int}, fbnd::Vector{Int}) where T
#   if isleaf(nd) 
#     _factor_leaf(A, nd.int, nd.bnd)
#   else
#     if !isnothing(nd.left)
#       left = _factor_branch(A, nd.left)
#       # figure out where the children indices go in the parent
#       left.fint = findfirst.(isequal.(nd.int), Ref(left.bnd))
#       left.fbnd = findfirst.(isequal.(nd.bnd), Ref(left.bnd))
#       intl = left.bnd[left.fint];
#       bndl = left.bnd[left.fbnd];
#     else
#       intl = Vector{Int}()
#       bndl = Vector{Int}()
#       left = nothing
#     end
#     if !isnothing(nd.right)
#       # modify computer father indices while going down
#       right = _factor_branch(A, nd.right)
#       right.fint = findfirst.(isequal.(nd.int), Ref(right.bnd))
#       right.fbnd = findfirst.(isequal.(nd.bnd), Ref(right.bnd))
#       intr = right.bnd[right.fint];
#       bndr = right.bnd[right.fbnd];
#     else
#       intl = Vector{Int}()
#       bndl = Vector{Int}()
#       right = nothing
#     end
#     # check that we really got all the degrees of freedom
#     int = [intl; intr]
#     bnd = [bndl; bndr]
#   end
# end

# function _factor_leaf(A::SparseMatrixCSC{T, Int}, int::Vector{Int}, bnd::Vector{Int}, fint::Vector{Int}, fbnd::Vector{Int}) where T
#   D = A[int, int]
#   L = Matrix(A[bnd, int]) / D # converts right/left-hand side to dense first
#   R = D \ Matrix(A[int, bnd])
#   S = A[bnd,bnd] - A[bnd, int] * R
#   SolverNode(D,S,L,R,int,bnd)
# end

# # creates the parent node but also updates the index sets of the children to indicate where boundary indices are found in the parent
# # maybe move this into the factorization process
# function parent!(left, right, inter, bound)
#   # compute where the children indices can be found in the parent
#   left.finter = findfirst.(isequal.(inter), Ref(left.bound))
#   left.fbound = findfirst.(isequal.(bound), Ref(left.bound))
#   right.finter = findfirst.(isequal.(inter), Ref(right.bound))
#   right.fbound = findfirst.(isequal.(bound), Ref(right.bound))
# end