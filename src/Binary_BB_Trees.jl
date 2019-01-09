__precompile__(true)

module Binary_BB_Trees

using Printf
using GeometryTypes: HomogenousMesh, Face, Point
using Rotations
using LinearAlgebra
using StaticArrays
using DataStructures
using Statistics
using NearestNeighbors
using NumericalTricks

include("box_types.jl")
include("mesh.jl")
include("util.jl")
include("vector_cache.jl")
include("tree_types.jl")
include("BB_intersection.jl")
include("blob_types.jl")
include("top_down.jl")
include("obb_construction.jl")


export
    # box_types.jl
    BoundingBox,
    AABB,
    OBB,
    boxArea,
    boxVolume,
    combine_BB,

    # mesh.jl
    Tri,
    Tet,
    eMesh,
    as_tet_eMesh,
    as_tri_eMesh,
    n_point,
    n_tri,
    n_tet,
    get_tri,
    get_tet,
    get_point,
    scale!,
    dh_transform_mesh!,
    mesh_inplace_rekey!,
    mesh_remove_unused_points!,
    delete_triangles!,
    mesh_repair!,
    output_eMesh_half_plane,
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
    get_vertices_32,
    get_faces_32,
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
    triTetMeshToTreeAABB,

    # top_down.jl

    # obb_construction.jl
    fit_tet_obb,
    fit_tri_obb,
    obb_tree_from_aabb

end # module
