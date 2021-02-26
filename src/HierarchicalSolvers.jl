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

  import Base.getproperty, Base.setproperty!, Base.size, Base.eltype, Base.getindex, Base.*, Base./, Base.\
  import HssMatrices.isleaf, HssMatrices.isbranch

  # HierarchicalSolvers.jl
  export SolverOptions
  # nesteddissection.jl
  export NDNode, NestedDissection, parse_elimtree, postorder, getinterior
  # blockmatrix.jl
  export BlockMatrix
  # factornode.jl
  export FactorNode, solve, solve!
  # factorization.jl
  export symfact!, factor

  mutable struct SolverOptions
    swlevel::Int
  end
  
  # set default values
  function SolverOptions(::Type{T}; args...) where T
    opts = SolverOptions(
      3               # swlevel
    )
    for (key, value) in args
      setfield!(opts, key, value)
    end
    opts
  end
  SolverOptions(; args...) = SolverOptions(Float64; args...)
  
  function copy(opts::SolverOptions; args...)
    opts_ = SolverOptions()
    for field in fieldnames(typeof(opts))
      setfield!(opts_, field, getfield(opts, field))
    end
    for (key, value) in args
      setfield!(opts_, key, value)
    end
    opts_
  end
  
  function chkopts!(opts::SolverOptions)
    opts.swlevel >= 0 || throw(ArgumentError("swlevel"))
  end

  include("nesteddissection.jl")
  include("blockmatrix.jl")
  include("factornode.jl")
  include("factorization.jl")
end # module
