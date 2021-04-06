### Definitions of datastructure and operators to store the nested dissection structure
# nested dissection is stored as a binary tree containing InvertedIndices
# The nested dissection datastructure contains not only indices to be eliminated but als the indives that they are connected to.
# Written by Boris Bonev, Feb. 2021

# Relies on the Binary tree structure defined in HssMatrices (Maybe change that?)
const UnsortedNestedDissection = BinaryNode{Tuple{Vector{Int}, Vector{Int}}}
const SortedNestedDissection = BinaryNode{Tuple{UnitRange{Int}, UnitRange{Int}}}

const NestedDissection = UnsortedNestedDissection

for nd_type in (:UnsortedNestedDissection,:SortedNestedDissection)
  @eval begin
    _getproperty(x::$nd_type, ::Val{s}) where {s} = getfield(x, s)
    _getproperty(x::$nd_type, ::Val{:int}) = x.data[1]
    _getproperty(x::$nd_type, ::Val{:bnd}) = x.data[2]
    getproperty(x::$nd_type, s::Symbol) = _getproperty(x, Val{s}())
    _setproperty!(x::$nd_type, ::Val{s}, a) where {s} = setfield!(x, s, a)
    _setproperty!(x::$nd_type, ::Val{:int}, a) = (x.data =  (a, x.data[2]))
    _setproperty!(x::$nd_type, ::Val{:bnd}, a) = (x.data =  (x.data[1], a))
    setproperty!(x::$nd_type, s::Symbol, a) = _setproperty!(x, Val{s}(), a)
  end
end

UNDNode(int::Vector{Int}, bnd::Vector{Int}) = BinaryNode((int, bnd))
UNDNode(int::Vector{Int}, bnd::Vector{Int}, left::UnsortedNestedDissection, right::UnsortedNestedDissection) = BinaryNode((int, bnd), left, right)
UNDNode(left::UnsortedNestedDissection, right::UnsortedNestedDissection) = BinaryNode((union(left.bnd, right.bnd), empty(left.bnd)), left, right)

## Symbolic factorization routine
# returns a reordered nessted dissection tree as well as a nested dissection tree containing the local indices
function symfact!(nd::UnsortedNestedDissection)
  _symfact!(nd, 1)
  nd_loc = symfact!(nd)
  nd_loc.int = collect(1:length(nd.bnd))
  nd_loc.bnd = Vector{Int}()
  return nd, nd_loc
end
function _symfact!(nd::UnsortedNestedDissection, level)
  if isleaf(nd)
    nd_loc = UNDNode(Vector{Int}(), Vector{Int}())
  else
    if !isnothing(nd.left)
      left_loc = _symfact!(nd.left, level+1)
      # figure out where the children indices go in the parent
      left_loc.int = findall(in(nd.int), nd.left.bnd)
      left_loc.bnd = findall(in(nd.bnd), nd.left.bnd)
      # reordered indices for the current node
      intl = nd.left.bnd[left_loc.int];
      bndl = nd.left.bnd[left_loc.bnd];
    else
      intl = Vector{Int}()
      bndl = Vector{Int}()
      left_loc = nothing
    end
    if !isnothing(nd.right)
      right_loc = _symfact!(nd.right, level+1)
      right_loc.int = findall(in(nd.int), nd.right.bnd)
      right_loc.bnd = findall(in(nd.bnd), nd.right.bnd)
      intr = nd.right.bnd[right_loc.int];
      bndr = nd.right.bnd[right_loc.bnd];
    else
      intr = Vector{Int}()
      bndr = Vector{Int}()
      right_loc = nothing
    end
    # check that we really got all the degrees of freedom
    nd.int = [intl; intr]
    nd.bnd = [bndl; bndr]
    nd_loc = UNDNode(Vector{Int}(), Vector{Int}(), left_loc, right_loc)
  end
  return nd_loc
end



## Other convenience functions
# get all indices in post-order
function postorder(nd::UnsortedNestedDissection)
  ind = Vector{Int}()
  for x in PostOrderDFS(nd)
    ind = [ind; x.int]
  end
  ind = [ind; nd.bnd]
end

# recursively compute the interior
getinterior(nd::UnsortedNestedDissection) = _getinterior!(nd, Vector{Int}())
function _getinterior!(nd::UnsortedNestedDissection, interior::Vector{Int})
  if !isnothing(nd.left) interior = _getinterior!(nd.left, interior) end
  if !isnothing(nd.right) interior = _getinterior!(nd.right, interior) end
  return [interior; nd.int]
end
getinterior(nd::SortedNestedDissection) = 1:nd.int[end]
getboundary(nd::NestedDissection) = nd.bnd

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
  snodes = Stack{UnsortedNestedDissection}()

  while !isempty(sind)
    i = first(sind)
    # moving up/down or a leaf?
    if rsons[i] == -1 && lsons[i] == -1 # at a leaf
      push!(snodes, UNDNode(inter[1:ninter[i], i], bound[1:nbound[i], i]))
      ilast = pop!(sind)
    elseif ilast == rsons[i] # moving up from the right
      right = pop!(snodes)
      if lsons[i] != -1
        left = pop!(snodes)
      else
        left = nothing
      end
      push!(snodes, UNDNode(inter[1:ninter[i], i], bound[1:nbound[i], i], left, right))
      ilast = pop!(sind)
    elseif ilast == lsons[i] && rsons[i] == -1 # moving up from the left but can't move down
      left = pop!(snodes)
      right = nothing
      push!(snodes, UNDNode(inter[1:ninter[i], i], bound[1:nbound[i], i], left, right))
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
