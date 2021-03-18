### 2x2 Block matrix structure. Convenience datastrucutre to hold our diagonal blocks for elimination.
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
# getindex(B::BlockMatrix, i::Int, j::Int) = getindex(B, [i], [j])[1]
getindex(B::BlockMatrix, i::Int, j::AbstractRange) = getindex(B, [i], j)[:]
getindex(B::BlockMatrix, i::AbstractRange, j::Int) = getindex(B, i, [j])[:]
getindex(B::BlockMatrix, i::AbstractRange, j::AbstractRange) = getindex(B, collect(i), collect(j))
getindex(B::BlockMatrix, ::Colon, ::Colon) = full(B)
getindex(B::BlockMatrix, i, ::Colon) = getindex(B, i, 1:size(B,2))
getindex(B::BlockMatrix, ::Colon, j) = getindex(B, 1:size(B,1), j)
function getindex(B::BlockMatrix, i::Vector{Int}, j::Vector{Int})
  m1, n1 = size(B.A11)
  A = Matrix{eltype(B)}(undef, length(i), length(j))
  if length(i) == 0 || length(j) == 0 return A end
  ii1 = i .<= m1; i1 = i[ii1]
  ii2 = i .> m1; i2 = i[ii2] .- m1
  jj1 = j .<= n1; j1 = j[jj1]
  jj2 = j .> n1; j2 = j[jj2] .- n1
  A[ii1, jj1] = B.A11[i1, j1]
  A[ii1, jj2] = B.A12[i1, j2]
  A[ii2, jj1] = B.A21[i2, j1]
  A[ii2, jj2] = B.A22[i2, j2]
  return A
end
function getindex(B::BlockMatrix, i::Int, j::Int)
  m1,n1 = size(B.A11)
  if !(1 ≤ i ≤ size(B,1)) || !(1 ≤ j ≤ size(B,2)); throw(BoundsError("Index out of bounds")); end
  if 1 ≤ i ≤ m1
    if 1 ≤ j ≤ n1
      return B.A11[i,j]
    else
      return B.A12[i,j-n1]
    end
  else
    if 1 ≤ j ≤ n1
      return B.A21[i-m1,j]
    else
      return B.A22[i-m1,j-n1]
    end
  end
end

# cast to matrix
Matrix(B::BlockMatrix) = [Matrix(B.A11) Matrix(B.A12); Matrix(B.A21) Matrix(B.A22)]

# adjoint/transpose
adjoint(B::BlockMatrix) = BlockMatrix(adjoint(B.A11), adjoint(B.A21), adjoint(B.A12), adjoint(B.A22))
transpose(B::BlockMatrix) = BlockMatrix(transpose(B.A11), transpose(B.A21), transpose(B.A12), transpose(B.A22))
copy(B::BlockMatrix) = BlockMatrix(copy(B.A11), copy(B.A12), copy(B.A21), copy(B.A22))

## Some basic linear algebra routines

# general fall-back routine
function *(A::BlockMatrix, B::AbstractMatrix)
  m1, n1 = size(A.A11)
  [A.A11*B[1:n1,:] .+ A.A12*B[n1+1:end,:]; A.A21*B[1:n1,:] .+ A.A22*B[n1+1:end,:]]
end
function *(A::AbstractMatrix, B::BlockMatrix)
  m1, n1 = size(B.A11)
  [A[:,1:m1]*B.A11 + A[:,m1+1:end]*B.A21 A[:,1:m1]*B.A12 + A[:,m1+1:end]*B.A22]
end
# specializations
function *(A::BlockMatrix, B::BlockMatrix)
  size(A,2) == size(B,1) ||  throw(DimensionMismatch("First dimension of B does not match second dimension of A. Expected $(size(A, 2)), got $(size(B, 1))"))
  size(A.A11,2) == size(B.A11,1) ||  throw(DimensionMismatch("Block rows of B do not match block columns of A. Expected $(size(A.A11, 2)), got $(size(B.A11, 1))"))
  BlockMatrix(A.A11*B.A11 + A.A12*B.A21, A.A11*B.A12 + A.A12*B.A22, A.A21*B.A11 + A.A22*B.A21, A.A21*B.A12 + A.A22*B.A22)
end
# specialize to make sure block matrices play nice with low-rank matrices
*(A::BlockMatrix, B::LowRankMatrix) = LowRankMatrix(A*B.U,B.V)
*(A::LowRankMatrix, B::BlockMatrix) = LowRankMatrix(A.U,(A.V'*B)')



# #\(A::BlockMatrix, B::AbstractMatrix) = blockldiv!(similar(B), A, copy(B))
# #/(A::AbstractMatrix, B::BlockMatrix) = blockrdiv!(similar(A), copy(A), B)

# specialized routines for computing A \ B overwriting B
blockldiv!(A::BlockMatrix, B::Matrix; atol::Float64, rtol::Float64) = B = blockldiv!(similar(B), A, B; atol, rtol)
function blockldiv!(Y::Matrix, A::BlockMatrix, B::Matrix; atol::Float64, rtol::Float64)
  m1,n1 = size(A.A11)
  if ishss(A.A11) && ishss(A.A12) && ishss(A.A21) && ishss(A.A22)
    S22 = A.A22 - A.A21*(A.A11\A.A12)
    S22 = recompress!(S22, atol=atol, rtol=rtol)
  elseif typeof(A.A11) <: HssMatrix && typeof(A.A22) <: HssMatrix
    error("Not implemented yet")
  else
    S22 = A.A22 .- A.A21*(A.A11\convert(Matrix, A.A12))
  end
  Y[1:n1,:] = A.A11\B[1:n1,:]
  Y[n1+1:end,:] = B[n1+1:end,:] - A.A21*Y[1:n1,:]
  Y[n1+1:end,:] = S22\Y[n1+1:end,:]
  Y[1:n1,:] = Y[1:n1,:] .- A.A11\(A.A12*Y[n1+1:end,:])
  return Y
end
# compute B / A overwriting B
blockrdiv!(A::Matrix, B::BlockMatrix; atol::Float64, rtol::Float64) = A = blockrdiv!(similar(A), A, B; atol, rtol)
function blockrdiv!(Y::Matrix, A::Matrix, B::BlockMatrix; atol::Float64, rtol::Float64)
  m1,n1 = size(B.A11)
  if ishss(B.A11) && ishss(B.A12) && ishss(B.A21) && ishss(B.A22)
    S22 = B.A22 - B.A21*(B.A11\B.A12)
    S22 = recompress!(S22, atol=atol, rtol=rtol)
  elseif ishss(B.A11) || ishss(B.A22)
    error("Not implemented yet")
  else
    S22 = B.A22 .- B.A21*(B.A11\convert(Matrix, B.A12))
  end
  #if ishss(B.A11); @infiltrate; end
  Y[:,1:m1] = A[:,1:m1]/B.A11
  Y[:,m1+1:end] = A[:,m1+1:end] - Y[:,1:m1]*B.A12
  Y[:,m1+1:end] = Y[:,m1+1:end]/S22
  Y[:,1:m1] = Y[:,1:m1] - (Y[:,m1+1:end]*B.A21)/B.A11
  return Y
end

function blockldiv(A::BlockMatrix, B::BlockMatrix)
  size(A,1) == size(B,1) || throw(DimensionMismatch("First dimension of A doesn't match first dimension of B. Expected $(size(B,1)), but got $(size(A,1))"))
  #size(A.A11,1) == size(B,1) || throw(DimensionMismatch("Block structure of A doesn't match first dimension of B. Expected $(size(B,1)), but got $(size(A,1))"))
  # Form the Schur complement
  if ishss(A.A11) && ishss(A.A12) && ishss(A.A21) && ishss(A.A22)
    S22 = A.A22 - A.A21*(A.A11\A.A12)
  elseif ishss(A.A11) || ishss(A.A22)
    error("Not implemented yet")
  else
    S22 = A.A22 - A.A21*(A.A11\convert(Matrix, A.A12))
  end
  # TODO: this can be optimized by storing factorizations
  # perform on the first column
  B11 = A.A11\B.A11
  B21 = B.A21 - A.A21*B11
  B21 = S22\B21
  B11 = B11 - A.A11\(A.A12*B21)
  # repeat for the second column
  B12 = A.A11\convert(Matrix, B.A12)
  B22 = B.A22 - A.A21*B12;
  B22 = S22\B22
  B12 = B12 - A.A11\(A.A12*B22); 
  return BlockMatrix(B11, B12, B21, B22)
end

function blockrdiv(B::BlockMatrix, A::BlockMatrix)
  # safety is off on this one
  # compute the Schur complement first
  if ishss(A.A11) && ishss(A.A12) && ishss(A.A21) && ishss(A.A22)
    S22 = A.A22 - A.A21*(A.A11\A.A12)
  elseif typeof(A.A11) <: HssMatrix && typeof(A.A2) <: HssMatrix
    error("Not implemented yet")
  else
    S22 = A.A22 - A.A21*(A.A11\convert(Matrix, A.A12))
  end
  # perform solve step on first block row
  B11 = B.A11/A.A11
  B12 = B.A12 - B11*A.A12
  B12 = B12/S22
  B11 = B11 - (B12*A.A21)/A.A11
  # repeat for second block row
  B21 = convert(Matrix, B.A21)/A.A11
  B22 = B.A22 - B21*A.A12
  B22 = B22/S22; 
  B21 = B21 - (B22*A.A21)/A.A11
  return BlockMatrix(B11, B12, B21, B22)
end