using LinearAlgebra, SparseArrays, LowRankApprox, Plots, HssMatrices
using BenchmarkTools

include("../src/HierarchicalSolvers.jl")
using .HierarchicalSolvers
import .HierarchicalSolvers: BlockMatrix

m1 = 120; n1 = 120;
m2 = 80; n2 = 80;
m = m1+m2; n = n1+n2

# create two matrices with identical blocking
A11 = randn(m1, n1); A12 = sprandn(m1, n2, 0.2); A21 = LowRankMatrix(randn(m2, 2), randn(n1, 2)); A22 = randn(m2,n2);
B11 = randn(m1, n1); B12 = LowRankMatrix(randn(m1, 10), randn(n2, 10)); B21 = sprandn(m2, n1, 0.1); B22 = randn(m2,n2);
A = BlockMatrix(A11, A12, A21, A22)
B = BlockMatrix(B11, B12, B21, B22)

println(norm(A*B - Matrix(A)*Matrix(B)))

A*sprandn(m,100,0.1)
B*sprandn(m,100,0.1)