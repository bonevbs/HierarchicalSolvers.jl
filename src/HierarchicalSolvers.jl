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

  import Base.getproperty

  const swlevel = 3

  # nesteddissection.jl
  export NDNode, NestedDissection, parse_elimtree, postorder

  include("nesteddissection.jl")
end # module
