### Definitions of the solver/preconditioner datastructure
# Written by Boris Bonev, Feb. 2021

mutable struct SolverNode{T<:Number} #<: Factorization
  D::Union{Matrix{T}, SparseMatrixCSC{T}}
  S::Union{Matrix{T}, HssMatrix{T}}
  L::Union{Matrix{T}, LowRankMatrix{T}}
  R::Union{Matrix{T}, LowRankMatrix{T}}

  # rename these to something more meaningful!!
  int::Vector{Int}
  bnd::Vector{Int}

  int_loc::Vector{Int}
  bnd_loc::Vector{Int}

  left::Union{SolverNode{T}, Nothing}
  right::Union{SolverNode{T}, Nothing}

  # internal constructors with checks for dimensions
  global function _SolverNode(D::Union{Matrix{T}, SparseMatrixCSC{T}}, S::Union{Matrix{T}, HssMatrix{T}},
    L::Union{Matrix{T}, LowRankMatrix{T}}, R::Union{Matrix{T}, LowRankMatrix{T}},
    int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int}) where T
    new{T}(D, S, L, R, int, bnd, int_loc, bnd_loc, nothing, nothing)
  end
  # parent constructor, finds the local indices of the children indices automatically 
  global function _SolverNode(D::Union{Matrix{T}, SparseMatrixCSC{T}}, S::Union{Matrix{T}, HssMatrix{T}},
    L::Union{Matrix{T}, LowRankMatrix{T}}, R::Union{Matrix{T}, LowRankMatrix{T}},
    int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int},
    left::SolverNode{T}, right::SolverNode{T}) where T
    # maybe also implement a check to make sure that disjointedness is guaranteed
    new{T}(D, S, L, R, int, bnd, int_loc, bnd_loc, left, right)
  end
end
SolverNode(D::Union{Matrix{T}, SparseMatrixCSC{T}}, S::Union{Matrix{T}, HssMatrix{T}},
L::Union{Matrix{T}, LowRankMatrix{T}}, R::Union{Matrix{T}, LowRankMatrix{T}},
int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int}) where T = _SolverNode(D, S, L, R, int, bnd, int_loc, bnd_loc)
SolverNode(D::Union{Matrix{T}, SparseMatrixCSC{T}}, S::Union{Matrix{T}, HssMatrix{T}},
L::Union{Matrix{T}, LowRankMatrix{T}}, R::Union{Matrix{T}, LowRankMatrix{T}},
int::Vector{Int}, bnd::Vector{Int}, int_loc::Vector{Int}, bnd_loc::Vector{Int},
left::SolverNode{T}, right::SolverNode{T}) where T = _SolverNode(D, S, L, R, int, bnd, int_loc, bnd_loc, left, right)

function _apply_factorization
end