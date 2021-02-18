### 2x2 Block matrix structure.
# Written by Boris Bonev, Feb. 2021

# 2x2 block matrix that can hold any type of 
struct BlockMatrix{T, T11 <: AbstractMatrix{T}, T12 <: AbstractMatrix{T}, T21 <: AbstractMatrix{T}, T22 <: AbstractMatrix{T}} <: AbstractMatrix{T}
  A11::T11
  A12::T12
  A21::T21
  A22::T22

  # inner constructor checks for consistency among dimensions
  function BlockMatrix(A11::AbstractMatrix{T}, A12::AbstractMatrix{T}, A21::AbstractMatrix{T}, A22::AbstractMatrix{T}) where T
    size(A11, 1) == size(A12, 1) || throw(DimensionMismatch("first dimension of A11 and A12 do not match. Expected $(size(A11, 1)), got $(size(A12, 1))"))
    size(A11, 2) == size(A21, 2) || throw(DimensionMismatch("second dimension of A11 and A12 do not match. Expected $(size(A11, 2)), got $(size(A21, 2))"))
    size(A22, 1) == size(A21, 1) || throw(DimensionMismatch("first dimension of A22 and A21 do not match. Expected $(size(A22, 1)), got $(size(A21, 1))"))
    size(A22, 2) == size(A12, 2) || throw(DimensionMismatch("second dimension of A22 and A12 do not match. Expected $(size(A22, 2)), got $(size(A12, 2))"))
    new{T, typeof(A11), typeof(A12), typeof(A21), typeof(A22)}(A11, A12, A21, A22)
  end
end

# basic overrides
eltype(::Type{BlockMatrix{T, S}}) where {T, S} = T
size(B::BlockMatrix) = size(B.A11) .+ size(B.A22)
size(B::BlockMatrix, dim::Int) = size(B)[dim]

# getindex
# this definition might be inefficient, but as it is mainly for display, this is fine
function _getindex(B::BlockMatrix, i::Int, j::Int)
  m1,n1 = size(B.A11)
  ret = zero(eltype(B))
  if 1 ≤ i ≤ m1
    if 1 ≤ j ≤ n1
      ret = B.A11[i,j]
    else
      ret = B.A12[i,j-n1]
    end
  else
    if 1 ≤ j ≤ n1
      ret = B.A12[i-m1,j]
    else
      ret = B.A22[i-m1,j-n1]
    end
  end
  return ret
end
function getindex(B::BlockMatrix, i::Int, j::Int)
  m,n = size(B)
  if 1 ≤ i ≤ m && 1 ≤ j ≤ n
      _getindex(B,i,j)
  else
      throw(BoundsError())
  end
end
getindex(B::BlockMatrix, i::Int, jr::AbstractRange) = transpose(eltype(B)[B[i,j] for j=jr])
getindex(B::BlockMatrix, ir::AbstractRange, j::Int) = eltype(B)[B[i,j] for i=ir]
getindex(B::BlockMatrix, ir::AbstractRange, jr::AbstractRange) = eltype(B)[B[i,j] for i=ir,j=jr]