__precompile__(true)

module Binary_BB_Trees

using GeometryTypes: HomogenousMesh, Face, Point
using Rotations
using LinearAlgebra
using StaticArrays
using DataStructures
using Statistics


include("util.jl")
include("box_types.jl")
include("vector_cache.jl")
include("tree_types.jl")
include("BB_intersection.jl")
# include("mesh_types.jl")
include("mesh_for_box.jl")
include("blob_types.jl")

export
    # util.jl
    svSvToAABB,
    extract_HomogenousMesh_face_vertices,

    # box_types.jl
    boundingBox,
    AABB,

    # vector_cache.jl
    VectorCache,
    expand!,
    returnNext,
    addCacheItem!,
    empty!,
    isempty,
    length,

    # tree_types.jl
    bin_BB_Tree,
    TT_Cache,
    update_TT_Cache!,
    length,
    is_leaf,
    is_not_leaf,
    treeDepth,
    leafNumber,
    tree_tree_intersect,
    extractData,

    # OBB,
    boxArea,
    boxVolume,
    boxMinMax,
    combineAABB,

    # intersection_tests.jl
    BB_BB_intersect,

    # # mesh_types.jl
    # SurfaceVolumeMesh,

    # mesh_for_box.jl
    outputBoxVolMesh,
    outputOrientedBoxFaces,
    basicBoxPoints,

    # blob_types.jl
    blob,
    triTetMeshToTreeAABB

end # module
