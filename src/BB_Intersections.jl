# we generate code in this module, so precompile where possible
__precompile__(true)

module BB_Intersections

using Rotations
using LinearAlgebra
using StaticArrays
using RigidBodyDynamics.Spatial

import Statistics

include("box_types.jl")
include("intersection_tests.jl")

export
    # box_types
    AABB,
    OBB,

    # intersection_tests
    BB_BB_intersect
end # module
