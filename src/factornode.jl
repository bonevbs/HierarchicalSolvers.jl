### Definitions of the solver/preconditioner datastructure
# Written by Boris Bonev, Feb. 2021

const MatFact{T} = Union{AbstractMatrix{T}, Factorization{T}}

# datastructure to hold the hiearchical factorization
mutable struct FactorNode{T<:Number, TD<:MatFact, TS<:AbstractMatrix{T}, TL<:AbstractMatrix{T}, TR<:AbstractMatrix{T}} <: Factorization{T}
  D::TD # diagonal block
  S::TS # Schur complement
  L::TL # left Gauss transform
  R::TR # right Gauss transform

  # rename these to something more meaningful!!
  int::Vector{Int}
  bnd::Vector{Int}

  int_loc::Vector{Int}
  bnd_loc::Vector{Int}

  # bindings to the children nodes
  left::Union{FactorNode{T}, Nothing}
  right::Union{FactorNode{T}, Nothing}

  # internal constructors with checks for dimensions
  global function _FactorNode(D::MatFact{T}, S::AbstractMatrix{T}, L::AbstractMatrix{T}, R::AbstractMatrix{T},
     int::AbstractVector{Int}, bnd::AbstractVector{Int}, int_loc::AbstractVector{Int}, bnd_loc::AbstractVector{Int}) where T
    new{T, typeof(D), typeof(S), typeof(L), typeof(R)}(D, S, L, R, int, bnd, int_loc, bnd_loc, nothing, nothing)
  end
  # parent constructor, finds the local indices of the children indices automatically 
  global function _FactorNode(D::MatFact{T}, S::AbstractMatrix{T}, L::AbstractMatrix{T}, R::AbstractMatrix{T},
      int::AbstractVector{Int}, bnd::AbstractVector{Int}, int_loc::AbstractVector{Int}, bnd_loc::AbstractVector{Int}, left::FactorNode{T}, right::FactorNode{T}) where T
    # maybe also implement a check to make sure that disjointedness is guaranteed
    new{T, typeof(D), typeof(S), typeof(L), typeof(R)}(D, S, L, R, int, bnd, int_loc, bnd_loc, left, right)
  end
end

# outer constructors
FactorNode(D::MatFact{T}, S::AbstractMatrix{T}, L::AbstractMatrix{T}, R::AbstractMatrix{T}, int, bnd, int_loc, bnd_loc) where T = _FactorNode(D, S, L, R, int, bnd, int_loc, bnd_loc)
FactorNode(D::MatFact{T}, S::AbstractMatrix{T}, L::AbstractMatrix{T}, R::AbstractMatrix{T}, int, bnd, int_loc, bnd_loc, left::FactorNode{T}, right::FactorNode{T}) where T = _FactorNode(D, S, L, R, int, bnd, int_loc, bnd_loc, left, right)

eltype(::FactorNode{T}) where T = T
Base.show(io::IO, node::FactorNode) = print(io, "FactorNode{$(eltype(node))}")

# usual checks to determine where we are
isleaf(F::FactorNode) = isnothing(F.left) && isnothing(F.right)
isbranch(F::FactorNode) = !isnothing(F.left) && !isnothing(F.right)

# function that recursively computes maximum rank
function maxrank(F::FactorNode)
  rkl = 0; rkr = 0
  if !isnothing(F.left) rkl = maxrank(F) end
  if !isnothing(F.right) rkr = maxrank(F) end
  rkS = ishss(F.S) ? hssrank(F.S) : 0
  rkL = F.L <: LowRankMatrix ? hssrank(F.L) : 0
  rkR = F.r <: LowRankMatrix ? hssrank(F.R) : 0
  return max(rkl, rkr, rkS, rkL, rkR)
end

## apply the factorization to a matrix
#solve(F, rhs) = solve!(F, copy(rhs))
#ldiv!(F::FactorNode, rhs::AbstractMatrix) = solve!(F, rhs)
ldiv!(F::FactorNode, B) = ldiv!(similar(B), F, B)
function ldiv!(C::AbstractVector{T}, F::FactorNode{T}, B::AbstractVector{T}) where T
  m = length(B)
  mC = ldiv!(reshape(C, m, 1), F, reshape(B, m, 1))
  C .= reshape(mC, m)
end
function ldiv!(C::AbstractMatrix{T}, F::FactorNode{T}, B::AbstractMatrix{T}) where T
  C .= B
  _lsolve!(F, C)
  _dsolve!(F, C)
  C[F.bnd, :] = F.S \ C[F.bnd, :]
  _rsolve!(F, C)
end

# recursive solution routines
function _lsolve!(F::FactorNode, rhs::AbstractMatrix)
  if !isnothing(F.left) rhs = _lsolve!(F.left, rhs) end
  if !isnothing(F.right) rhs = _lsolve!(F.right, rhs) end
  rhs[F.bnd,:] = rhs[F.bnd,:] - F.L*rhs[F.int,:]
  return rhs
end
function _rsolve!(F::FactorNode, rhs::AbstractMatrix)
  rhs[F.int,:] = rhs[F.int,:] - F.R*rhs[F.bnd,:]
  if !isnothing(F.left) rhs = _rsolve!(F.left, rhs) end
  if !isnothing(F.right) rhs = _rsolve!(F.right, rhs) end
  return rhs
end
function _dsolve!(F::FactorNode, rhs::AbstractMatrix)
  if !isnothing(F.left) rhs = _dsolve!(F.left, rhs) end
  if !isnothing(F.right) rhs = _dsolve!(F.right, rhs) end
  #rhs[F.int,:] .= F.D\rhs[F.int,:]
  if isa(F.D, BlockFactorization)
    rhs[F.int,:] = blockldiv!(F.D, rhs[F.int,:])
  else
    rhs[F.int,:] = F.D\rhs[F.int,:]
  end
  return rhs
end

# # functions for nice access
# function \(F::FactorNode{TF}, B::AbstractVector{TB}) where {TF,TB}
#   T = promote_type(TF, TB)
#   FT = convert(FactorNode{T}, F)
#   BT = (T == TB ? B : convert(Array{T}, B))
#   CT = Array{T}(undef, size(A,2))
#   ldiv!(CT, AT, BT)
# end
# function \(F::FactorNode{TF}, B::AbstractMatrix{TB}) where {TF,TB}
#   T = promote_type(TA, TB)
#   FT = convert(FactorNode{T}, F)
#   BT = (T == TB ? B : convert(Array{T}, B))
#   CT = Array{T}(undef, size(A,2), size(B,2))
#   ldiv!(CT, AT, BT)
# end