
findMinSVSV(a::SVector{3,Float64}, b::SVector{3,Float64}) = min.(a, b)
findMinSVSV(sv::SVector{3, SVector{3,Float64}}) = min.(min.(sv[1], sv[2]),      sv[3])
findMinSVSV(sv::SVector{4, SVector{3,Float64}}) = min.(min.(sv[1], sv[2]), min.(sv[3], sv[4]))

findMaxSVSV(a::SVector{3,Float64}, b::SVector{3,Float64}) = max.(a, b)
findMaxSVSV(sv::SVector{3, SVector{3,Float64}}) = max.(max.(sv[1], sv[2]),      sv[3])
findMaxSVSV(sv::SVector{4, SVector{3,Float64}}) = max.(max.(sv[1], sv[2]), max.(sv[3], sv[4]))

function minMaxToCenterExtent(box_min::SVector{3,Float64}, box_max::SVector{3,Float64})
    center = (box_max + box_min) * 0.5
    extent = (box_max - box_min) * 0.5
    return center, extent
end

function svSvToAABB(vert::SVector{N, SVector{3,Float64}}) where {N}
    box_min = findMinSVSV(vert)
    box_max = findMaxSVSV(vert)
    center, extent = minMaxToCenterExtent(box_min, box_max)
    return AABB(center, extent)
end

function aabb_for_points(a::SVector{3,Float64}, b::SVector{3,Float64})
    min_ = findMinSVSV(a, b)
    max_ = findMaxSVSV(a, b)
    center, extent = minMaxToCenterExtent(min_, max_)
    return AABB(center, extent)
end

function combine_BB(a::BB_Type, b::BB_Type) where {BB_Type <: BoundingBox}
    # TODO: put in utility

    min_1, max_1 = calc_min_max(a)
    min_2, max_2 = calc_min_max(b)
    aabb = svSvToAABB(SVector{4,SVector{3,Float64}}(min_1, max_1, min_2, max_2))
    return BB_Type(aabb)
end

calc_min_max(a::AABB) = a.c - a.e, a.c + a.e

function calc_min_max(point::Vector{SVector{3,T}}) where {T}
    min_val = SVector{3,Float64}(+Inf, +Inf, +Inf)
    max_val = SVector{3,Float64}(-Inf, -Inf, -Inf)
    for k = 1:length(point)
        point_k = point[k]
        min_val = min.(min_val, point_k)
        max_val = max.(max_val, point_k)
    end
    return min_val, max_val
end

function calc_min_max(hm::HomogenousMesh)
    point, _ = extract_HomogenousMesh_face_vertices(hm)
    return calc_min_max(point)
end

# function calc_min_max(sv::SVector{N,Float64}) where {N}
function calc_min_max(sv::SVector{3,Float64})
    min_ = findMinSVSV(sv)
    max_ = findMaxSVSV(sv)
    return min_, max_
end

function calc_min_max(a::SVector{3,Float64}, b::SVector{3,Float64})
    min_ = findMinSVSV(a, b)
    max_ = findMaxSVSV(a, b)
    return min_, max_
end

function calc_aabb(arg_in)
    min_val, max_val = calc_min_max(arg_in)
    center, extent = minMaxToCenterExtent(min_val, max_val)
    return AABB(center, extent)
end

function calc_aabb(a::SVector{3,Float64}, b::SVector{3,Float64})
    min_val, max_val = calc_min_max(a, b)
    center, extent = minMaxToCenterExtent(min_val, max_val)
    return AABB(center, extent)
end

function sortEdgeFace(v::SVector{3,Int64}, k::Int64)
    three = 3
    (1 <= k <= three) || error("a triangle has three sides")
    i1 = v[mod1(k + 1, three)]
    i2 = v[mod1(k + 2, three)]
    return SVector{2,Int64}(minmax(i1, i2))
end

function sortEdgeFace(v::SVector{4,Int64}, k::Int64)
    four = 4
    (1 <= k <= four) || error("a tetrahedron has four sides")
    i1 = v[mod1(k + 1, four)]
    i2 = v[mod1(k + 2, four)]
    i3 = v[mod1(k + 3, four)]
    i1, i2 = minmax(i1, i2)  # --> i2 is not lowest
    i2, i3 = minmax(i2, i3)  # --> i3 is highest
    i1, i2 = minmax(i1, i2)  # --> i1 and i2 are sorted
    return SVector{3,Int64}(i1, i2, i3)
end

function extract_HomogenousMesh_face_vertices(hm::HomogenousMesh)
    point = [SVector{3,Float64}(k) for k = hm.vertices]
    vec_tri = [SVector{3,Int64}(k) for k = hm.faces]
    return point, vec_tri
end

get_h_mesh_vertices(hm::HomogenousMesh) = [SVector{3,Float64}(k) for k = hm.vertices]
get_h_mesh_faces(hm::HomogenousMesh) = [SVector{3,Int64}(k) for k = hm.faces]

get_h_mesh_vertices_32(hm::HomogenousMesh) = [Point{3,Float32}(k) for k = hm.vertices]
get_h_mesh_faces_32(hm::HomogenousMesh) = [Face{3,Int32}(k) for k = hm.faces]

get_vertices_32(e_mesh::eMesh{Tri,T2}) where {T2} = [Point{3,Float32}(k) for k = e_mesh.point]
get_faces_32(e_mesh::eMesh{Tri,T2}) where {T2} = [Face{3,Int32}(k) for k = e_mesh.tri]

# function scale_HomogenousMesh!(mesh::HomogenousMesh, scale::Float64)
#     return scale_HomogenousMesh!(mesh, scale * ones(SVector{3,Float64}))
# end
# function scale_HomogenousMesh!(mesh::HomogenousMesh, scale::SVector{3,Float64})
#     for k = 1:length(mesh.vertices)
#         mesh.vertices[k] = mesh.vertices[k] .* scale
#     end
#     return nothing
# end
#
# function transform_HomogenousMesh!(mesh::HomogenousMesh; rot::Rotation=one(Quat{Float64}), trans::SVector{3,Float64}=zeros(SVector{3,Float64}))
#     R = RotMatrix(rot)
#     for k = 1:length(mesh.vertices)
#         mesh.vertices[k] = R * mesh.vertices[k] + trans
#     end
#     return nothing
# end
#
# function repair_mesh(mesh_leaky::HomogenousMesh)
#     # Currently this function: merges duplicate points in defective meshes.
#     # Deletes unused points
#
#     function shortest_side(abc::SVector{3,SVector{3,Float64}})
#         a, b, c = abc
#         return min(norm(a - b), norm(b - c), norm(c - a))
#     end
#
#     # All points closer to each other than hald a side length are duplicates
#     ind = get_h_mesh_faces(mesh_leaky)
#     point = get_h_mesh_vertices(mesh_leaky)
#     balltree = BallTree(point; reorder = false)
#     idxs, dists = knn(balltree, point, 20, true)
#     min_side_length = minimum([shortest_side(point[k]) for k = ind])
#     idxs = inrange(balltree, point, min_side_length * 0.499)
#
#     # Determine mapping of old indices to new indices
#     n_point = length(point)
#     i_orig_prime = collect(1:n_point)
#     for k = 1:n_point
#         k_current = i_orig_prime[k]
#         i_orig_prime[idxs[k]] .= k_current
#     end
#     unique_i_orig_prime = unique(i_orig_prime)
#     i_prime_new = zeros(Int64, n_point) .- 9999
#     i_prime_new[unique_i_orig_prime] .= collect(1:length(unique_i_orig_prime))
#
#     # Create repaired mesh
#     point_new = point[unique_i_orig_prime]
#     i_orig_new = i_prime_new[i_orig_prime]
#     ind_new = [i_orig_new[k] for k = ind]
#     return HomogenousMesh(faces=ind_new, vertices=point_new)
# end

function recursivly_rotate!(i::Vector{SVector{3,Int64}}, p::Vector{SVector{3,Float64}}, i_detele::BitArray{1},
        k::Int64, n̂::SVector{3,Float64}, d::Float64, perm::Int64, is_hard::Bool)

    function util(ind::Int64)
        point = p[ind]
        proj_val = dot(point, n̂) + d
        return point, proj_val, proj_val >= 0.0
    end

    (perm == 4) && error("perm is 4")
    i1 = i[k][perm]
    i2 = i[k][mod1(perm + 1, 3)]
    i3 = i[k][mod1(perm + 2, 3)]
    p1, v1, b1 = util(i1)
    p2, v2, b2 = util(i2)
    p3, v3, b3 = util(i3)
    if b2 != b3
        recursivly_rotate!(i, p, i_detele, k, n̂, d, perm + 1, is_hard)
    else
        le_sum = b1 + b2 + b3
        if le_sum == 3
        elseif is_hard || (le_sum == 0)
            i_detele[k] = true
        else
            c_12 = weightPoly(p1, p2, v1, v2)
            c_31 = weightPoly(p3, p1, v3, v1)
            push!(p, c_12)
            i_12 = length(p)
            push!(p, c_31)
            i_31 = length(p)
            if le_sum == 1
                i[k] = SVector{3,Int64}(i1, i_12, i_31)
            elseif le_sum == 2
                i[k] = SVector{3,Int64}(i2, i3, i_12)
                push!(i, SVector{3,Int64}(i3, i_31, i_12))
            end
        end
    end
end

# function crop_mesh(mesh::HomogenousMesh, n̂::SVector{3,Float64}, d::Float64, is_hard::Bool=false)
#     m = deepcopy(mesh)
#     p = get_h_mesh_vertices(m)
#     i = get_h_mesh_faces(m)
#     n_faces_orig = length(i)
#     bool_delete = falses(n_faces_orig)
#     area_ = [area(p[k]) for k = i]
#     min_area = minimum(area_)
#     min_area_tol = min_area / 100
#     for k = 1:n_faces_orig
#         recursivly_rotate!(i, p, bool_delete, k, n̂, d, 1, is_hard)
#     end
#     area_ = [area(p[k]) for k = i]
#     i_small = findall(area_ .<= min_area_tol)
#     i_delete = findall(bool_delete)
#     append!(i_delete, i_small)
#     i_delete = sort(unique(i_delete))
#
#     deleteat!(i, i_delete)
#     new_hm = HomogenousMesh(vertices=p, faces=i)
#     return repair_mesh(new_hm)
# end
