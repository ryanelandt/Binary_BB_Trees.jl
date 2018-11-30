__precompile__(true)

module Binary_BB_Trees

using GeometryTypes: HomogenousMesh, Face, Point
using Rotations
using LinearAlgebra
using StaticArrays
using DataStructures
using Statistics
using NearestNeighbors
using NumericalTricks

include("box_types.jl")
include("util.jl")
include("vector_cache.jl")
include("tree_types.jl")
include("BB_intersection.jl")
include("blob_types.jl")
include("top_down.jl")
include("mesh.jl")


export
    # box_types.jl
    AABB,
    boxArea,
    boxVolume,
    combineAABB,

    # mesh.jl
    Tri,
    Tet,
    eMesh,
    as_tet_eMesh,
    as_tri_eMesh,
    n_points,
    n_tri,
    n_tet,
    mesh_inplace_rekey!,
    mesh_remove_unused_points!,
    delete_triangles!,
    mesh_repair!,
    output_eMesh_box,
    output_eMesh_hole,

    # util.jl
    svSvToAABB,
    calc_min_max,
    calc_aabb,
    extract_HomogenousMesh_face_vertices,
    get_h_mesh_vertices,
    get_h_mesh_faces,
    get_h_mesh_vertices_32,
    get_h_mesh_faces_32,
    scale_HomogenousMesh!,
    transform_HomogenousMesh!,
    repair_mesh,
    crop_mesh,

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

    # intersection_tests.jl
    BB_BB_intersect,

    # blob_types.jl
    blob,
    triTetMeshToTreeAABB

    # top_down.jl

end # module
