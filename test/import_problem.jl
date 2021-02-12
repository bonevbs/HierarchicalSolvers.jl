include("../src/HierarchicalSolvers.jl")
using .HierarchicalSolvers
using HssMatrices
using MAT
using LinearAlgebra, SparseArrays, Plots
using DataStructures

## Read the problem from MAT file
file = matopen("./test/test.mat")
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
inter   = convert(Matrix{Int}, elim_tree["ninter"]);
bound   = convert(Matrix{Int}, elim_tree["nbound"]);
# just convert everything into the appropriate formats


parse_nested_dissection(fathers, lsons, rsons, ninter, inter, nbound, bound)

# create the nested dissection from the elimination tree DataStructures
#nested_dissection(fathers, lsons, rsons, ninter, inter, nbound, bound)

# function first_integer(fathers::Vector{Int})
#   println(a[1])
# end