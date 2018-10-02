# Binary_BB_Trees

[![Build Status](https://travis-ci.com/ryanelandt/Binary_BB_Trees.jl.svg?branch=master)](https://travis-ci.com/ryanelandt/Binary_BB_Trees.jl)
[![codecov.io](https://codecov.io/github/ryanelandt/Binary_BB_Trees.jl/coverage.svg?branch=master)](https://codecov.io/github/ryanelandt/Binary_BB_Trees.jl?branch=master)

This module implements BVH (bounding volume hierarchy) creation and intersection tests for binary trees made up of AABBs (axis aligned bounding boxes).



### Create AABB Tree

Load a mesh using [MeshIO](https://github.com/JuliaIO/MeshIO.jl).
Assuming that the loaded mesh is called `mesh` one creates an AABB tree by calling:
```Julia
tree_mesh = triTetMeshToTreeAABB(mesh)
```

### Determine Intersecting AABB/AABB Pairs

```Julia
using LinearAlgebra
using StaticArrays

tt_cache = TT_Cache()
update_TT_Cache!(tt_cache, SVector{3,Float64}(0,0,0), SMatrix{3,3,Float64,9}(I))
tree_tree_intersect(tt_cache, tree_mesh, tree_mesh)
```
