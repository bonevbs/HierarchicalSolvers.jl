using LinearAlgebra, SparseArrays, LowRankApprox, Plots
using BenchmarkTools

include("../src/HierarchicalSolvers.jl")
using .HierarchicalSolvers
using HssMatrices

m1 = 120; n1 = 120;
m2 = 80; n2 = 80;

# create two matrices with identical blocking
A11 = randn(m1, n1); A12 = sprandn(m1, n2, 0.2); A21 = LowRankMatrix(randn(m2, 20), randn(n1, 20)); A22 = randn(m2,n2);
B11 = randn(m1, n1); B12 = LowRankMatrix(randn(m1, 10), randn(n2, 10)); B21 = sprandn(m2, n1, 0.1); B22 = randn(m2,n2);
A = BlockMatrix(A11, A12, A21, A22)
B = BlockMatrix(B11, B12, B21, B22)