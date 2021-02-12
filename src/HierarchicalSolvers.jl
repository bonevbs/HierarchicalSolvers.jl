### HierarchicalSolvers.jl module
# A simple Julia package implementing general NestedDissection type solvers and our hierarchical preconditioner
# Written by Boris Bonev, Feb. 2021
__precompile__()
module HierarchicalSolvers

  using AbstractTrees
  using HssMatrices
  using LinearAlgebra
  using SparseArrays
  using DataStructures

  const swlevel = 3

  # nesteddissection.jl
  export NDNode, NestedDissection, parse_nested_dissection

  include("nesteddissection.jl")
end # module
