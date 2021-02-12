### Definitions of datastructures and basic constructors and operators
# Written by Boris Bonev, Feb. 2021

# nested dissection is stored as a binary tree containing InvertedIndices
const NestedDissection = BinaryNode{Tuple{Vector{Int}, Vector{Int}}}

function getindex(A::NestedDissection, d::Symbol)
  if     d == :int  return A.data[1]
  elseif d == :bnd  return A.data[2]
  else            throw(KeyError(d))
  end
end

NDNode(int::Vector{Int}, bnd::Vector{Int}) = BinaryNode((int, bnd))
NDNode(int::Vector{Int}, bnd::Vector{Int}, left::NestedDissection, right::NestedDissection) = BinaryNode((int, bnd), left, right)
NDNode(left::NestedDissection, right::NestedDissection) = BinaryNode((union(left.bnd, right.bnd), empty(left.bnd)), left, right)


# function postorder(nd::NestedDissection)
#   if !isnothing(node.left)
# end

# function for reading in 
function parse_nested_dissection(fathers::Vector{Int}, lsons::Vector{Int}, rsons::Vector{Int}, ninter::Vector{Int}, inter::Matrix{Int}, nbound::Vector{Int}, bound::Matrix{Int})
  nnodes = length(fathers)
  nnodes == length(lsons) == length(rsons) == length(niter) == length(nbound) == size(inter,2) == size(bound,2) || throw(DimensionMismatch("dimensions inconsistent among inputs"))

  # find the root
  iroot = findall(x->x==-1, fathers);
  if length(iroot) != 1 throw(ArgumentError("found either less than or more than one root.")) end
  i = iroot[1]
  s = Stack{Int}()
  snd = Stack{NDNode}()

  while !isempty(s)
    # moving up/down or a leaf?
    if rsons[i] == -1 && lsons[i] == -1
      push!(snd, NDNode(inter[1:ninter[i], i], bound[1:nbound[i], i]))
      i = pop!(s)
    # moving up from the left nad a right child remains
    elseif ilast == lsons[i] && rsons[i] != -1
      ndleft = ndlast
      ilast = 1
      i = pop!(s)
    elseif ilast == rsons[i] || (ilast == lsons[i] && rsons[i] == -1)
      ndlast = NDNode(inter[1:ninter[i], i], bound[1:nbound[i], i], ndleft, ndright)
      ilast = i
    else
      push!(s, rsons[i])
      push!(s, lsons[i])
      ilast = i
      i = lsons[i]
    end
  end
end
