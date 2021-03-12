### Some extensions to LowRankMatrix.. this should be moved to LowRankApprox.jl
# Written by Boris Bonev, Feb. 2021

function mul!(C::AbstractMatrix, A::LowRankApprox.LowRankMatrix, B::AbstractMatrix, α::Number, β::Number)
  size(A,2) == size(B,1) ||  throw(DimensionMismatch("First dimension of B does not match second dimension of A. Expected $(size(A, 2)), got $(size(B, 1))"))
  size(C) == (size(A,1), size(B,2)) ||  throw(DimensionMismatch("Dimensions of C don't match up with A and B."))
  tmp = mul!(similar(B,rank(A),size(B,2)), A.V', B)
  return mul!(C, A.U, tmp, α, β)
end