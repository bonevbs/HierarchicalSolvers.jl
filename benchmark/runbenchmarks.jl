### import_problem.jl
# Example file to load a sample problem from a .mat-file, where the nested dissection structure is described in a special, serialized formats
#
# Written by Boris Bonev, Feb. 2021

using LinearAlgebra, SparseArrays, IterativeSolvers, Plots
using MAT

#using HssMatrices
include("../src/HierarchicalSolvers.jl")
using .HierarchicalSolvers

using BenchmarkTools

## Read the problem from MAT file
file = matopen("./benchmark/benchmark_problem.mat")
elim_tree = read(file, "elim_tree") # note that this does NOT introduce a variable ``varname`` into scope
A = read(file, "A")
b = read(file, "b")
close(file)

# enforce tree datastructure to be in the form of Int vectors
fathers = convert(Vector{Int}, dropdims(elim_tree["fathers"], dims=1));
lsons   = convert(Vector{Int}, dropdims(elim_tree["lsons"],   dims=1));
rsons   = convert(Vector{Int}, dropdims(elim_tree["rsons"],   dims=1));
ninter  = convert(Vector{Int}, dropdims(elim_tree["ninter"],  dims=1));
nbound  = convert(Vector{Int}, dropdims(elim_tree["nbound"],  dims=1));
inter   = convert(Matrix{Int}, elim_tree["inter"]);
bound   = convert(Matrix{Int}, elim_tree["bound"]);
# just convert everything into the appropriate formats

nd = parse_elimtree(fathers, lsons, rsons, ninter, inter, nbound, bound)

F = factor(A, nd, swlevel = 0, atol=1e-3, rtol=1e-3);
x = solve!(F, copy(b));

println("Benchmarking hierarchical factorization...")
@btime F = factor(A, nd, swlevel = 0, atol=1e-3, rtol=1e-3)

println("Benchmarking solve step...")
@btime x = solve!(F, copy(b))
