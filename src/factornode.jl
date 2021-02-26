### Definitions of the solver/preconditioner datastructure
# Written by Boris Bonev, Feb. 2021

# datastructure to hold the hiearchical factorization
mutable struct FactorNode{T<:Number} <: Factorization{T}
  D::Union{Matrix{T}, SparseMatrixCSC{T}, BlockMatrix{T}}
  S::Union{Matrix{T}, HssMatrix{T}}
  L::Union{Matrix{T}, LowRankMatrix{T}}
  R::Union{Matrix{T}, LowRankMatrix{T}}

  # rename these to something more meaningful!!
  int::Vector{Int}
  bnd::Vector{Int}

  int_loc::Vector{Int}
  bnd_loc::Vector{Int}

  left::Union{FactorNode{T}, Nothing}
  right::Union{FactorNode{T}, Nothing}

  # internal constructors with checks for dimensions
  global function _FactorNode(D::Union{Matrix{T}, SparseMatrixCSC{T}, BlockMatrix{T}}, S::Union{Matrix{T}, HssMatrix{T}},
    L::Union{Matrix{T}, LowRankMatrix{T}}, R::Union{Matrix{T}, LowRankMatrix{T}},
    int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int}) where T
    new{T}(D, S, L, R, int, bnd, int_loc, bnd_loc, nothing, nothing)
  end
  # parent constructor, finds the local indices of the children indices automatically 
  global function _FactorNode(D::Union{Matrix{T}, SparseMatrixCSC{T}, BlockMatrix{T}}, S::Union{Matrix{T}, HssMatrix{T}},
    L::Union{Matrix{T}, LowRankMatrix{T}}, R::Union{Matrix{T}, LowRankMatrix{T}},
    int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int},
    left::FactorNode{T}, right::FactorNode{T}) where T
    # maybe also implement a check to make sure that disjointedness is guaranteed
    new{T}(D, S, L, R, int, bnd, int_loc, bnd_loc, left, right)
  end
end

# outer constructors
FactorNode(D::Union{Matrix{T},SparseMatrixCSC{T}, BlockMatrix{T}}, S::Union{Matrix{T}, HssMatrix{T}},
L::Union{Matrix{T}, LowRankMatrix{T}}, R::Union{Matrix{T}, LowRankMatrix{T}},
int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int}) where T = _FactorNode(D, S, L, R, int, bnd, int_loc, bnd_loc)
FactorNode(D::Union{Matrix{T},SparseMatrixCSC{T}, BlockMatrix{T}}, S::Union{Matrix{T}, HssMatrix{T}},
L::Union{Matrix{T}, LowRankMatrix{T}}, R::Union{Matrix{T}, LowRankMatrix{T}},
int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int},
left::FactorNode{T}, right::FactorNode{T}) where T = _FactorNode(D, S, L, R, int, bnd, int_loc, bnd_loc, left, right)

eltype(::FactorNode{T}) where T = T

# usual checks to determine where we are
isleaf(F::FactorNode) = isnothing(F.left) && isnothing(F.right)
isbranch(F::FactorNode) = !isnothing(F.left) && !isnothing(F.right)

## apply the factorization to a matrix
solve(F, rhs) = solve!(F, copy(rhs))
function solve!(F::FactorNode, rhs::AbstractMatrix)
  rhs = _lsolve!(F, rhs)
  rhs = _dsolve!(F, rhs)
  rhs[F.bnd, :] .= F.S \ rhs[F.bnd, :]
  rhs = _rsolve!(F, rhs)
end

# recursive solution routines
function _lsolve!(F::FactorNode, rhs::AbstractMatrix)
  if !isnothing(F.left) rhs = _lsolve!(F.left, rhs) end
  if !isnothing(F.right) rhs = _lsolve!(F.right, rhs) end
  rhs[F.bnd,:] .= rhs[F.bnd,:] .- F.L*rhs[F.int,:]
  return rhs
end
function _rsolve!(F::FactorNode, rhs::AbstractMatrix)
  rhs[F.int,:] .= rhs[F.int,:] .- F.R*rhs[F.bnd,:]
  if !isnothing(F.left) rhs = _rsolve!(F.left, rhs) end
  if !isnothing(F.right) rhs = _rsolve!(F.right, rhs) end
  return rhs
end
function _dsolve!(F::FactorNode, rhs::AbstractMatrix)
  if !isnothing(F.left) rhs = _dsolve!(F.left, rhs) end
  if !isnothing(F.right) rhs = _dsolve!(F.right, rhs) end
  rhs[F.int,:] .= F.D\rhs[F.int,:]
  return rhs
end