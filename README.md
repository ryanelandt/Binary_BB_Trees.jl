# Binary_BB_Trees

[![Build Status](https://travis-ci.com/ryanelandt/Binary_BB_Trees.jl.svg?branch=master)](https://travis-ci.com/ryanelandt/Binary_BB_Trees.jl)
[![codecov.io](https://codecov.io/github/ryanelandt/Binary_BB_Trees.jl/coverage.svg?branch=master)](https://codecov.io/github/ryanelandt/Binary_BB_Trees.jl?branch=master)

This module implements BVH (bounding volume hierarchy) creation and intersection tests for binary trees made up of AABBs (axis aligned bounding boxes).
This package was designed to perform binary collision detection for dynamics software.
On a high level this package does two things:
1. create AABB trees for triangular and tetrahedral meshes and
2. determine all intersections between two AABB trees for an arbitrary relative position and orientation.

### Create AABB Tree

Load a mesh using [MeshIO](https://github.com/JuliaIO/MeshIO.jl).
Assuming that the loaded mesh is called `mesh` one creates an AABB tree by calling:
```Julia
tree_mesh = triTetMeshToTreeAABB(mesh)
```

Tetrahedral meshes can be created by passing a vector of points and vector of indices to `triTetMeshToTreeAABB` (see below).
```Julia
function triTetMeshToTreeAABB(point::Vector{SVector{3,Float64}}, vec_tri_tet::Vector{SVector{N,Int64}}) where {N}
```

### Determine Intersecting AABB/AABB Pairs

Intersecting pairs are tested using the OBB/OBB SAT test described in Real-Time Collision Detection by Ericson.

```Julia
using LinearAlgebra
using StaticArrays

tt_cache = TT_Cache()  # the indices of intersecting pairs are put in here
relative_translation = SVector{3,Float64}(0,0,0)
relative_orientation = SMatrix{3,3,Float64,9}(I)
update_TT_Cache!(tt_cache, relative_translation, relative_orientation)  # update relative position and orientation
tree_tree_intersect(tt_cache, tree_mesh, tree_mesh)
```

### Tree Construction Algorithm

An intersection query needs to return all intersecting pairs.
A good AABB tree accomplishes this goal while testing few pairs.
Heuristics that govern BHV creation generally try to do some combination of the following:
1. balance tree,
2. minimize surface area and
3. minimize volume.

The bottom up algorithm implemented in this module uses a heuristic that considers these three things.
The heuristic in invariant to translation and scaling by a constant.

The code bottom up code in this module is somewhat obtuse, but necessary to construct trees in `n log(n)` time.
This paragraph attempts to explain how it works in a way that people can understand.
The algorithm starts with all bounding boxes separated.
The number of leaves in each box, denoted by `n`, is 1 when the algorithm starts.
The cost of each box is tabulated according to a heuristic that assigns cost of the form `a * volume + b * surface_area + n * log2(2*n)`; the purpose for the `n * log2(2*n)` will be made clear later soon.
Each box corresponds to a triangle/tetrahedron.
The algorithm calculates the marginal cost of combining boxes that belong to neighbors (i.e. triangles that touch each other).
If the `n * log2(2*n)` term were not present, there would be no penalty for creating an unbalanced tree.
The algorithm merges neighbors one at a time in the manner that increases the sum of all costs the least until either a single box remains, or no possible merges are left (mesh is not connected).
