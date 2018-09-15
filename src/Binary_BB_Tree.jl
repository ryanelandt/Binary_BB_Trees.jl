__precompile__(true)

module Binary_BB_Tree

using Rotations
using LinearAlgebra
using StaticArrays
using RigidBodyDynamics.Spatial

import Statistics

include("box_types.jl")
include("intersection_tests.jl")
include("tree_types.jl")

export
    # box_types
    AABB,
    OBB,
    bondingBox,

    # intersection_tests
    BB_BB_intersect

end # module
