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
include("blob_types.jl")
include("top_down.jl")


export
    # util.jl
    svSvToAABB,
    calc_min_max,
    calc_aabb,
    extract_HomogenousMesh_face_vertices,
    find_vector_point_AABB,
    get_h_mesh_vertices,
    get_h_mesh_faces,
    get_h_mesh_vertices_32,
    get_h_mesh_faces_32,
    scale_HomogenousMesh!,
    transform_HomogenousMesh!,

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

    # blob_types.jl
    blob,
    triTetMeshToTreeAABB

    # top_down.jl

end # module
