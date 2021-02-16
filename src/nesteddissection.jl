### Definitions of datastructure and operators to store the nested dissection structure
# nested dissection is stored as a binary tree containing InvertedIndices
# The nested dissection datastructure contains not only indices to be eliminated but als the indives that they are connected to.
# Written by Boris Bonev, Feb. 2021

# Relies on the Binary tree structure defined in HssMatrices (Maybe change that?)
const NestedDissection = BinaryNode{Tuple{Vector{Int}, Vector{Int}}}

_getproperty(x::NestedDissection, ::Val{s}) where {s} = getfield(x, s)
_getproperty(x::NestedDissection, ::Val{:int}) = x.data[1]
_getproperty(x::NestedDissection, ::Val{:bnd}) = x.data[2]
Base.getproperty(x::NestedDissection, s::Symbol) = _getproperty(x, Val{s}())
_setproperty!(x::NestedDissection, ::Val{s}, a) where {s} = setfield!(x, s, a)
_setproperty!(x::NestedDissection, ::Val{:int}, a) = (x.data =  (a, x.data[2]))
_setproperty!(x::NestedDissection, ::Val{:bnd}, a) = (x.data =  (x.data[1], a))
Base.setproperty!(x::NestedDissection, s::Symbol, a) = _setproperty!(x, Val{s}(), a)

NDNode(int::Vector{Int}, bnd::Vector{Int}) = BinaryNode((int, bnd))
NDNode(int::Vector{Int}, bnd::Vector{Int}, left::NestedDissection, right::NestedDissection) = BinaryNode((int, bnd), left, right)
NDNode(left::NestedDissection, right::NestedDissection) = BinaryNode((union(left.bnd, right.bnd), empty(left.bnd)), left, right)

# get all indices in post-order
function postorder(nd::NestedDissection)
  ind = Vector{Int}()
  for x in PostOrderDFS(nd)
    ind = [ind; x.int]
  end
  ind = [ind; nd.bnd]
end

# convenience function for reading in my own serialized elimination tree format
# probably of little to no use to others
function parse_elimtree(fathers::Vector{Int}, lsons::Vector{Int}, rsons::Vector{Int}, ninter::Vector{Int}, inter::Matrix{Int}, nbound::Vector{Int}, bound::Matrix{Int})
  nnodes = length(fathers)
  nnodes == length(lsons) == length(rsons) == length(ninter) == length(nbound) == size(inter,2) == size(bound,2) || throw(DimensionMismatch("dimensions inconsistent among inputs"))

  # find the root
  roots = findall(x->x==-1, fathers);
  if length(roots) != 1 throw(ArgumentError("found either less than or more than one root.")) end

  # prepare Stack for tree indices and for the assembly
  sind = Stack{Int}()
  push!(sind, roots[1])
  ilast = -2
  snodes = Stack{NestedDissection}()

  while !isempty(sind)
    i = first(sind)
    # moving up/down or a leaf?
    if rsons[i] == -1 && lsons[i] == -1 # at a leaf
      push!(snodes, NDNode(inter[1:ninter[i], i], bound[1:nbound[i], i]))
      ilast = pop!(sind)
    elseif ilast == rsons[i] # moving up from the right
      right = pop!(snodes)
      if lsons[i] != -1
        left = pop!(snodes)
      else
        left = nothing
      end
      push!(snodes, NDNode(inter[1:ninter[i], i], bound[1:nbound[i], i], left, right))
      ilast = pop!(sind)
    elseif ilast == lsons[i] && rsons[i] == -1 # moving up from the left but can't move down
      left = pop!(snodes)
      right = nothing
      push!(snodes, NDNode(inter[1:ninter[i], i], bound[1:nbound[i], i], left, right))
      ilast = pop!(sind)
    elseif (ilast == lsons[i] && rsons[i] != -1) || (lsons[i] == -1 && rsons != -1) # move down to the right
      ilast = i
      push!(sind, rsons[i])
    else # go down to the left
      ilast = i
      push!(sind, lsons[i])
    end
  end
  pop!(snodes)
end
