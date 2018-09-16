__precompile__(true)

module Binary_BB_Tree

using Rotations
using LinearAlgebra
using StaticArrays
using RigidBodyDynamics.Spatial

import Statistics

include("box_types.jl")
include("BB_intersection.jl")
include("vector_cache.jl")
include("tree_types.jl")
include("mesh_types.jl")
include("mesh_for_box.jl")

export
    # box_types.jl
    AABB,
    OBB,
    boundingBox,

    # intersection_tests.jl
    BB_BB_intersect,

    # vector_cache.jl
    vectorCache,
    expand!,
    returnNext,
    addCacheItem!,
    empty!,
    isempty,
    length,

    # tree_types.jl
    bin_BB_Tree,
    TT_Cache,
    isleaf,
    tree_tree_intersect,

    # mesh_types.jl
    SurfaceVolumeMesh,

    # mesh_for_box.jl
    boxMesh

end # module
