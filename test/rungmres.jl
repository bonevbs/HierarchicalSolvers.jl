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
#A, b, nd = read_problem("./test/poisson2d_p1_h64.mat")
A, b, nd = read_problem("./test/poisson3d_p1_h16_nmax10.mat")
nd, nd_loc = symfact!(nd)
perm = postorder(nd)
A = permute(A, perm, perm)
nd = permuted!(nd, invperm(perm))

pdeg=1
bsz = 60*(pdeg+1)*(pdeg+2)*(pdeg+3)

println("Problem parameters:")
println("   $(size(A)) matrix")
println("   $(depth(nd))-level nested-dissection")

println("Replacing sparse getindex")
@eval SparseArrays include("../src/mygetindex.jl");

println("Computing factorization without compression...")
Fa, ta = @timed factor(A, nd, nd_loc; swlevel = 0)
println("took ", ta, " seconds")
#@time Fa = factor(A, nd; swlevel = 0)
#xa = ldiv!(Fa, copy(b));
#println("rel. error without compression ", norm(A*xa-b)/norm(A\b))

println("Computing factorization with compression...")
Fc, tc = @timed factor(A, nd, nd_loc; swlevel = -2, swsize=4*bsz, atol=1e-3, rtol=1e-3, kest=200, stepsize=100, leafsize=bsz, verbose=true)
println("took ", tc, " seconds")
#println("Computing approximate solution...")
#xc = ldiv!(Fc, copy(b));
#@time xc = ldiv!(Fc, copy(b));
#println("rel. error with compression ", norm(A*xc-b)/norm(A\b))

# compare with fill-in reduction
x1, ch1 = gmres(A, b; Pr=Fa, reltol=1e-9, restart=30, log=true, maxiter=30)
x2, ch2 = gmres(A, b; Pr=Fc, reltol=1e-9, restart=30, log=true, maxiter=30)

plot(yaxis=:log)
plot!(ch1[:resnorm], marker=true, label="direct solver")
plot!(ch2[:resnorm], marker=true, label="hierarchical preconditioner")