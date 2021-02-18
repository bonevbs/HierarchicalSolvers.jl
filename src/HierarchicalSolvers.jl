### HierarchicalSolvers.jl module
# A simple Julia package implementing general NestedDissection type solvers and our hierarchical preconditioner
# Written by Boris Bonev, Feb. 2021
__precompile__()
module HierarchicalSolvers

  using AbstractTrees
  using HssMatrices
  using LinearAlgebra
  using SparseArrays
  using LowRankApprox
  using DataStructures

  import Base.getproperty
  import HssMatrices.isleaf, HssMatrices.isbranch

  const swlevel = 3

  # nesteddissection.jl
  export NDNode, NestedDissection, parse_elimtree, postorder
  # factornode.jl
  export FactorNode, solve, solve!
  # factorization.jl
  export symfact!, factor

  include("nesteddissection.jl")
  include("factornode.jl")
  include("factorization.jl")
end # module
