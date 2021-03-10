# utility routine to read in problems in the specific Matlab format that I use
# Written by Boris Bonev, March 2021
using MAT

function read_problem(filepath)
  ## Read the problem from MAT file
  file = matopen(filepath)
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

  # conver to nested dissection
  nd = parse_elimtree(fathers, lsons, rsons, ninter, inter, nbound, bound)
  return A, b, nd
end