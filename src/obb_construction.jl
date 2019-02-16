
function tet_perm_by_num(n::Int64)
    (n == 1) && (return SVector{4,Int64}(2,4,3,1))
    (n == 2) && (return SVector{4,Int64}(4,1,3,2))
    (n == 3) && (return SVector{4,Int64}(1,4,2,3))
    (n == 4) && (return SVector{4,Int64}(1,2,3,4))
    error("n must be 1:4")
end

get_projections(p::SVector{3,SVector{3,Float64}}, ê) = SVector(dot(p[1], ê), dot(p[2], ê), dot(p[3], ê))
get_projections(p::SVector{4,SVector{3,Float64}}, ê) = SVector(dot(p[1], ê), dot(p[2], ê), dot(p[3], ê), dot(p[4], ê))

function make_obb(p::SVector{N,SVector{3,Float64}}, i_start::Int64) where {N}
    (N == 3) || (N == 4) || error("must have 3 or 4 points")
    ê_1 = normalize(p[mod1(i_start + 1, 3)] - p[i_start])
    ê_3 = triangleNormal(p[1], p[2], p[3])
    ê_2 = cross(ê_3, ê_1)
    proj_1 = get_projections(p, ê_1)
    proj_2 = get_projections(p, ê_2)
    proj_3 = get_projections(p, ê_3)
    p_min = SVector(minimum(proj_1), minimum(proj_2), minimum(proj_3))
    p_max = SVector(maximum(proj_1), maximum(proj_2), maximum(proj_3))
    c, e = minMaxToCenterExtent(p_min, p_max)
    R = hcat(ê_1, ê_2, ê_3)
    return OBB(R * c, e, R)
end

fit_tri_obb(p::SVector{3,SVector{3,Float64}}) = make_obb(p, 1)  # all bounding boxes have the same surface area
function fit_tet_obb(p::SVector{4,SVector{3,Float64}}, ϵ_tet::SVector{4,Float64})
    (0.0 < volume(p)) || error("inverted tet")
    i_core = findfirst(ϵ_tet .!= 0.0)  # TODO: handle this more optimally
    p = p[tet_perm_by_num(i_core)]  # rearrange so that a plane is made with the first three
    obb_1 = make_obb(p, 1)
    obb_2 = make_obb(p, 2)
    obb_3 = make_obb(p, 3)
    area_1 = boxArea(obb_1)
    area_2 = boxArea(obb_2)
    area_3 = boxArea(obb_3)
    (max(        area_2, area_3) <= area_1) && (return obb_1)
    (max(area_1,         area_3) <= area_2) && (return obb_2)
    (max(area_1, area_2        ) <= area_3) && (return obb_3)
end

function obb_tree_from_aabb(e_tri::bin_BB_Tree, all_obb_tri::Vector{OBB})
    if is_leaf(e_tri)
        id_ = e_tri.id
        return bin_BB_Tree(id_, all_obb_tri[id_])
    else
        box = OBB(e_tri.box)
        node_1 = obb_tree_from_aabb(e_tri.node_1, all_obb_tri)
        node_2 = obb_tree_from_aabb(e_tri.node_2, all_obb_tri)
        return bin_BB_Tree(box, node_1, node_2)
    end
end










# function rearrange_longest_side(p1::SVector{3,Float64}, p2::SVector{3,Float64}, p3::SVector{3,Float64})
#     function tri_find_longest_side(p1::SVector{3,Float64}, p2::SVector{3,Float64}, p3::SVector{3,Float64})
#         d_32 = norm(p3 - p2)
#         d_13 = norm(p1 - p3)
#         d_21 = norm(p2 - p1)
#         return findmax([d_32, d_13, d_21])
#     end
#
#     for _ = 1:2
#         length_longest, i_side = tri_find_longest_side(p1, p2, p3)
#         if i_side != 3
#             p1, p2, p3 = p2, p3, p1
#         end
#     end
#     return p1, p2, p3
# end
#
# function create_triangle_R(p1::SVector{3,Float64}, p2::SVector{3,Float64}, p3::SVector{3,Float64})
#     dir_3 = triangleNormal(p1, p2, p3)
#     dir_1 = normalize(p2 - p1)
#     dir_2 = cross(dir_1, -dir_3)
#     return hcat(dir_1, dir_2, dir_3)
# end
#
# function tet_extent(p1::SVector{3,Float64}, p2::SVector{3,Float64}, p3::SVector{3,Float64}, R::SMatrix{3,3,Float64,9})
#     return tet_extent(p1, p2, p3, p1, R)
# end
# function tet_extent(p1::SVector{3,Float64}, p2::SVector{3,Float64}, p3::SVector{3,Float64}, p4::SVector{3,Float64},
#         R::SMatrix{3,3,Float64,9})
#
#     extent_1 = dot(p2 - p1, R[:,1]) / 2
#     extent_2 = dot(p3 - p1, R[:,2]) / 2
#     extent_3 = dot(p4 - p1, R[:,3]) / 2
#     extent = SVector{3,Float64}(extent_1, extent_2, extent_3)
#     @assert(all(-1.0e-15 <= extent[1]), "extent_1 is negative")
#     @assert(all(-1.0e-15 <= extent[2]), "extent_3 is negative")
#     @assert(all(-1.0e-15 <= extent[3]), "extent_3 is negative")
#     return extent
# end
#
# function fit_tri_obb(p_tri::SVector{3,SVector{3,Float64}})
#     p1, p2, p3 = p_tri
#
#     p1, p2, p3 = rearrange_longest_side(p1, p2, p3)
#     R = create_triangle_R(p1, p2, p3)
#     extent = tet_extent(p1, p2, p3, R)
#     center_w = R * extent + p1
#
#     return OBB(center_w, extent, R)
# end
#
# function tet_perm_by_num(n::Int64)
#     if n == 1
#         return SVector{4,Int64}(2,4,3,1)
#     elseif n == 2
#         return SVector{4,Int64}(4,1,3,2)
#     elseif n == 3
#         return SVector{4,Int64}(1,4,2,3)
#     elseif n == 4
#         return SVector{4,Int64}(1,2,3,4)
#     else
#         error("n must be 1:4")
#     end
# end
#
# function project_tet_top_into_plane(p1, p2, p3, p4)
#     # plane described by n_hat * x + d = 0
#     n_hat = triangleNormal(p1, p2, p3)
#     d = -dot(n_hat, p1)
#     the_proj = dot(n_hat, p4) + d
#     return p4 - the_proj * n_hat
# end
#
# function shoelace(x_::SVector{N,Float64}, y_::SVector{N,Float64}) where {N}
#     function shoelace(x__::SVector{N,Float64}, y__::SVector{N,Float64}, k::Int64) where {N}
#         k_plus = mod1(k + 1, N)
#         return x__[k] * y__[k_plus] - x__[k_plus] * y__[k]
#     end
#
#     cum_area = 0.0
#     for k = 1:N
#         cum_area += shoelace(x_, y_, k)
#     end
#     return cum_area
# end
#
# function four_point_convex_hull(p_proj::SVector{4,SVector{3,Float64}})
#     (sum(abs.([k[3] - p_proj[1][3] for k = p_proj])) < 1.0e-14) || error("last element is not similar")
#
#     # using Combinatorics
#     # the_combo = collect(permutations([1,2,3,4]))
#     # for (count, k) = enumerate(the_combo)
#     #     println("SVector{4,Int64}(", k[1], ", ", k[2], ", ", k[3], ", ", k[4], ifelse(count == 24, ")", "),"))
#     # end
#     # for (count, k) = enumerate(the_combo)
#     #     println("SVector{3,Int64}(", k[1], ", ", k[2], ", ", k[3], ifelse(count == 24, ")", "),"))
#     # end
#     combo_4 = (
#         SVector{4,Int64}(1, 2, 3, 4),
#         SVector{4,Int64}(1, 2, 4, 3),
#         SVector{4,Int64}(1, 3, 2, 4),
#         SVector{4,Int64}(1, 3, 4, 2),
#         SVector{4,Int64}(1, 4, 2, 3),
#         SVector{4,Int64}(1, 4, 3, 2),
#         SVector{4,Int64}(2, 1, 3, 4),
#         SVector{4,Int64}(2, 1, 4, 3),
#         SVector{4,Int64}(2, 3, 1, 4),
#         SVector{4,Int64}(2, 3, 4, 1),
#         SVector{4,Int64}(2, 4, 1, 3),
#         SVector{4,Int64}(2, 4, 3, 1),
#         SVector{4,Int64}(3, 1, 2, 4),
#         SVector{4,Int64}(3, 1, 4, 2),
#         SVector{4,Int64}(3, 2, 1, 4),
#         SVector{4,Int64}(3, 2, 4, 1),
#         SVector{4,Int64}(3, 4, 1, 2),
#         SVector{4,Int64}(3, 4, 2, 1),
#         SVector{4,Int64}(4, 1, 2, 3),
#         SVector{4,Int64}(4, 1, 3, 2),
#         SVector{4,Int64}(4, 2, 1, 3),
#         SVector{4,Int64}(4, 2, 3, 1),
#         SVector{4,Int64}(4, 3, 1, 2),
#         SVector{4,Int64}(4, 3, 2, 1)
#     )
#     combo_3 = (
#         SVector{3,Int64}(1, 2, 3),
#         SVector{3,Int64}(1, 2, 4),
#         SVector{3,Int64}(1, 3, 2),
#         SVector{3,Int64}(1, 3, 4),
#         SVector{3,Int64}(1, 4, 2),
#         SVector{3,Int64}(1, 4, 3),
#         SVector{3,Int64}(2, 1, 3),
#         SVector{3,Int64}(2, 1, 4),
#         SVector{3,Int64}(2, 3, 1),
#         SVector{3,Int64}(2, 3, 4),
#         SVector{3,Int64}(2, 4, 1),
#         SVector{3,Int64}(2, 4, 3),
#         SVector{3,Int64}(3, 1, 2),
#         SVector{3,Int64}(3, 1, 4),
#         SVector{3,Int64}(3, 2, 1),
#         SVector{3,Int64}(3, 2, 4),
#         SVector{3,Int64}(3, 4, 1),
#         SVector{3,Int64}(3, 4, 2),
#         SVector{3,Int64}(4, 1, 2),
#         SVector{3,Int64}(4, 1, 3),
#         SVector{3,Int64}(4, 2, 1),
#         SVector{3,Int64}(4, 2, 3),
#         SVector{3,Int64}(4, 3, 1),
#         SVector{3,Int64}(4, 3, 2)
#     )
#
#     x = SVector{4,Float64}([k[1] for k = p_proj])
#     y = SVector{4,Float64}([k[2] for k = p_proj])
#     area_3 = [shoelace(x[k], y[k]) for k = combo_3]
#     area_4 = [shoelace(x[k], y[k]) for k = combo_4]
#     val_3, ind_3 = findmax(area_3)
#     val_4, ind_4 = findmax(area_4)
#     return (val_4 <= val_3) ? combo_3[ind_3] : combo_4[ind_4]
# end
#
# function hull_points(i_hull::SVector{N,Int64}, k::Int64) where {N}
#     return i_hull[k], i_hull[mod1(k + 1, N)]
# end
#
# function find_bb_area(p_proj, n_hat, i_hull::SVector{N,Int64}, k::Int64) where {N}
#     R = create_R(p_proj, n_hat, i_hull, k)
#     aabb = calc_aabb(Vector([R' * k for k = p_proj]))
#     return aabb.e[1] * aabb.e[2]
# end
#
# function create_R(p_four::SVector{4,SVector{3,Float64}}, n_hat::SVector{3,Float64}, i_hull::SVector{N,Int64}, k::Int64) where {N}
#     i1, i2 = hull_points(i_hull, k)
#     dir_1 = normalize(p_four[i2] - p_four[i1])
#     dir_2 = cross(n_hat, dir_1)
#     return hcat(dir_1, dir_2, n_hat)
# end
#
# function bound_tet_by_first_three_plane(p_tet::SVector{4,SVector{3,Float64}}, i_core::Int64)
#     offset_o = centroid(p_tet)
#     p_tet = SVector{4,SVector{3,Float64}}([k - offset_o for k = p_tet])
#
#     p1_o, p2_o, p3_o, p4_o = p_tet[Binary_BB_Trees.tet_perm_by_num(i_core)]  # rotate until the core is upward
#     (0.0 < volume(p1_o, p2_o, p3_o, p4_o)) || error("negative volume")
#
#     n_hat = triangleNormal(p1_o, p2_o, p3_o)
#     p1_o, p2_o, p3_o = Binary_BB_Trees.rearrange_longest_side(p1_o, p2_o, p3_o)
#     p4_proj = project_tet_top_into_plane(p1_o, p2_o, p3_o, p4_o)
#
#     R_arbitrary = Binary_BB_Trees.create_triangle_R(p1_o, p2_o, p3_o)
#     p_plane_o = SVector{4,SVector{3,Float64}}(p1_o, p2_o, p3_o, p4_proj)
#     p_plane_n = SVector{4,SVector{3,Float64}}([R_arbitrary' * k for k = p_plane_o])
#     i_hull = four_point_convex_hull(p_plane_n)
#
#     ### FIND optimal bounding box for the inplane hull
#     areas = [find_bb_area(p_plane_n, n_hat, i_hull, k) for k = 1:length(i_hull)]
#     v_min, i_min = findmin(areas)
#     i1, i2 = hull_points(i_hull, i_min)
#
#     ### CREATE rotation matrix
#     R = create_R(p_plane_o, n_hat, i_hull, i_min)
#     p_o = SVector{4,SVector{3,Float64}}(p1_o, p2_o, p3_o, p4_o)
#     aabb = calc_aabb(Vector([R' * k for k = p_o]))
#
#     center_ = R * aabb.c + offset_o
#     return OBB(center_, aabb.e, R)
# end
#
# function fit_tet_obb(p_tet::SVector{4,SVector{3,Float64}}, ϵ_tet::SVector{4,Float64})
#     i_core = findfirst(ϵ_tet .!= 0.0)
#     return bound_tet_by_first_three_plane(p_tet, i_core)
# end
