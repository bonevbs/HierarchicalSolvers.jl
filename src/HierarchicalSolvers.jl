### HierarchicalSolvers.jl module
# A simple Julia package implementing general NestedDissection type solvers and our hierarchical preconditioner
# Written by Boris Bonev, Feb. 2021
__precompile__()
module HierarchicalSolvers

  using AbstractTrees
  using LinearAlgebra
  using HssMatrices
  using SparseArrays
  using LowRankApprox
  using DataStructures

  import Base: getproperty, setproperty!, size, eltype, getindex, *, /, \, copy, adjoint, transpose
  import LinearAlgebra: ldiv!, rdiv!, mul!
  import HssMatrices: isleaf, isbranch

  # HierarchicalSolvers.jl
  export SolverOptions
  # nesteddissection.jl
  export NDNode, NestedDissection, parse_elimtree, postorder, getinterior, getboundary
  # lowrankextensions.jl
  # blockmatrix.jl
  # factornode.jl
  export FactorNode, solve, solve!
  # factorization.jl
  export symfact!, factor

  mutable struct SolverOptions
    swlevel::Int
    swsize::Int
    atol::Float64
    rtol::Float64
    c_tol::Float64
    leafsize::Int
    kest::Int
    stepsize::Int
    verbose::Bool
  end
  
  # set default values
  function SolverOptions(; args...)
    opts = SolverOptions(
      5,      # switching level at which to start compression
      1000,   # minimum size for compression ## not implemented yet
      1e-6,   # absolute compression tolerance
      1e-6,   # relative compression tolerance
      0.5,    # relative factor to tune low-rank compression tolerance w.r.t HSS compression tolerance
      32,     # HSS leaf size
      30,     # rank estimate
      10,     # stepsize
      false,  # output Important information
      )
    for (key, value) in args
      setfield!(opts, key, value)
    end
    return opts
  end
  #SolverOptions(; args...) = SolverOptions(Float64; args...)
  
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
    opts.swsize ≥ 1 || throw(ArgumentError("swsize"))
    opts.atol ≥ 0. || throw(ArgumentError("atol"))
    opts.rtol ≥ 0. || throw(ArgumentError("rtol"))
    0. < opts.c_tol ≤ 1. || throw(ArgumentError("c_tol"))
    opts.leafsize ≥ 1 || throw(ArgumentError("leafsize"))
  end

  include("nesteddissection.jl")
  include("lowrankextensions.jl")
  include("blockmatrix.jl")
  include("factornode.jl")
  include("factorization.jl")
end # module
