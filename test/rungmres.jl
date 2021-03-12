### import_problem.jl
# Example file to load a sample problem from a .mat-file, where the nested dissection structure is described in a special, serialized formats
#
# Written by Boris Bonev, Feb. 2021

using LinearAlgebra, SparseArrays, IterativeSolvers, Plots, Random
Random.seed!(123);

include("../src/HierarchicalSolvers.jl")
using .HierarchicalSolvers
using HssMatrices

include("../util/read_problem.jl")
#A, b, nd = read_problem("./test/test.mat")
A, b, nd = read_problem("./test/poisson2d_p1_h64.mat")
println("Read in $(size(A)) matrix.")

println("Computing factorization without compression...")
Fa = factor(A, nd, swlevel = 0)
@time Fa = factor(A, nd, swlevel = 0)
xa = ldiv!(Fa, copy(b));
println("rel. error without compression ", norm(A*xa-b)/norm(A\b))

println("Computing factorization with compression...")
Fc = factor(A, nd, swlevel = -2, atol=1e-6, rtol=1e-6)
@profview Fc = factor(A, nd, swlevel = -2, atol=1e-6, rtol=1e-6)
xc = ldiv!(Fc, copy(b));
println("rel. error with compression ", norm(A*xc-b)/norm(A\b))

# pind = postorder(nd)
# spy(A[pind, pind])

# compare with fill-in reduction
x1, ch1 = gmres(A, b; Pr=Fa, reltol=1e-9, restart=20, log=true, maxiter=5)
x2, ch2 = gmres(A, b; Pr=Fc, reltol=1e-9, restart=20, log=true, maxiter=5)

plot(yaxis=:log)
plot!(ch1, :resnorm, marker=true, label="direct solver")
plot!(ch2, :resnorm, marker=true, label="hierarchical preconditioner")