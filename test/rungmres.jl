### import_problem.jl
# Example file to load a sample problem from a .mat-file, where the nested dissection structure is described in a special, serialized formats
#
# Written by Boris Bonev, Feb. 2021

using LinearAlgebra, SparseArrays, IterativeSolvers, Plots, Random
Random.seed!(123);

include("../src/HierarchicalSolvers.jl")
using .HierarchicalSolvers
using HssMatrices

#using Profile
#using Traceur
#using ProfileView
#using Cthulhu
using TimerOutputs
#using BenchmarkTools

include("../util/read_problem.jl")
#A, b, nd = read_problem("./test/test.mat")
#A, b, nd = read_problem("./test/poisson2d_p1_h64.mat")
A, b, nd = read_problem("./test/poisson2d_p1_h512.mat")
nd, nd_loc = symfact!(nd)
perm = postorder(nd)
A = permute(A, perm, perm)
nd = permuted!(nd, invperm(perm))

bsz = 60

println("Problem parameters:")
println("   $(size(A)) matrix")
println("   $(depth(nd))-level nested-dissection")

reset_timer!(HierarchicalSolvers.to)
println("Computing factorization without compression...")
Fa = factor(A, nd, nd_loc; swlevel = 0)
show(HierarchicalSolvers.to)
println()
@time factor(A, nd, nd_loc; swlevel = 0)
#@time Fa = factor(A, nd; swlevel = 0)
#xa = ldiv!(Fa, copy(b));
#println("rel. error without compression ", norm(A*xa-b)/norm(A\b))

reset_timer!(HierarchicalSolvers.to)
println("Computing factorization with compression...")
Fc = factor(A, nd, nd_loc; swlevel = -6, swsize=4*bsz, atol=1e-6, rtol=1e-6, kest=40, stepsize=10, leafsize=bsz, verbose=false)
show(HierarchicalSolvers.to)
println()
@time factor(A, nd, nd_loc; swlevel = -6, swsize=4*bsz, atol=1e-6, rtol=1e-6, kest=40, stepsize=10, leafsize=bsz, verbose=false)
#@trace( Fc = factor(A, nd; swlevel = -4, atol=1e-6, rtol=1e-6, kest=40, stepsize=10, leafsize=bsz, verbose=false), maxdepth=2)
#@profview Fc = factor(A, nd; swlevel = -4, atol=1e-6, rtol=1e-6, kest=40, stepsize=10, leafsize=bsz, verbose=false)
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