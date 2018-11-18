findMinSVSV(sv::SVector{3, SVector{3,Float64}}) = min.(min.(sv[1], sv[2]),      sv[3])
findMinSVSV(sv::SVector{4, SVector{3,Float64}}) = min.(min.(sv[1], sv[2]), min.(sv[3], sv[4]))
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

function calc_min_max(sv::SVector{N,Float64}) where {N}
    min_ = findMinSVSV(sv)
    max_ = findMaxSVSV(sv)
    return min_, max_
end

function calc_aabb(arg_in)
    min_val, max_val = calc_min_max(arg_in)
    center, extent = minMaxToCenterExtent(min_val, max_val)
    return AABB(center, extent)
end

function find_vector_point_AABB(point::Vector{SVector{3,Float64}})
    # min_val = point[1] * 1
    # max_val = point[2] * 2
    # for k = 1:length(point)
    #     point_k = point[k]
    #     min_val = min.(min_val, point_k)
    #     max_val = max.(max_val, point_k)
    # end
    min_val, max_val = calc_min_max(point)
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

function scale_HomogenousMesh!(mesh::HomogenousMesh, scale::Float64)
    return scale_HomogenousMesh!(mesh, scale * ones(SVector{3,Float64}))
end
function scale_HomogenousMesh!(mesh::HomogenousMesh, scale::SVector{3,Float64})
    for k = 1:length(mesh.vertices)
        mesh.vertices[k] = mesh.vertices[k] .* scale
    end
    return nothing
end

function transform_HomogenousMesh!(mesh::HomogenousMesh; rot::Rotation=one(Quat{Float64}), trans::SVector{3,Float64}=zeros(SVector{3,Float64}))
    R = RotMatrix(rot)
    for k = 1:length(mesh.vertices)
        mesh.vertices[k] = R * mesh.vertices[k] + trans
    end
    return nothing
end
