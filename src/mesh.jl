
abstract type Tri end
abstract type Tet end

struct eMesh{T1<:Union{Nothing,Tri},T2<:Union{Nothing,Tet}}
    point::Vector{SVector{3,Float64}}
    tri::Union{Nothing,Vector{SVector{3,Int64}}}
    tet::Union{Nothing,Vector{SVector{4,Int64}}}
    ϵ::Union{Nothing,Vector{Float64}}
    function eMesh( point::Vector{SVector{3,Float64}},
                    tri::Union{Nothing,Vector{SVector{3,Int64}}},
                    tet::Union{Nothing,Vector{SVector{4,Int64}}}=nothing,
                    ϵ::Union{Nothing,Vector{Float64}}=nothing)

        T1_ = ifelse(tri == nothing, Nothing, Tri)
        T2_ = ifelse(tet == nothing, Nothing, Tet)
        if T2_ == Tet
            @assert(isa(ϵ, Vector{Float64}))
            @assert(length(ϵ) == length(point), "length(ϵ) = $(length(ϵ)) but length(point) = $(length(point))")
            (length(ϵ) != 0) && @assert(maximum(ϵ) == 0.0, "strain must be zero on the surface of the volume mesh")
            for k = 1:length(tet)
                (0.0 < volume(point[tet[k]])) || error("inverted tetrahedron")
            end
        else
            @assert(ϵ == nothing)
        end
        (T1_ == T2_ == Nothing) && error("a whole lot of nothing")
        return new{T1_,T2_}(point, tri, tet, ϵ)
    end
    function eMesh(hm::HomogenousMesh, tet::Union{Nothing,Vector{SVector{4,Int64}}}=nothing,
            ϵ::Union{Nothing,Vector{Float64}}=nothing)

        point = [SVector{3,Float64}(k) for k = hm.vertices]
        tri = [SVector{3,Int64}(k) for k = hm.faces]
        return eMesh(point, tri, tet, ϵ)
    end
    eMesh{Tri,Nothing}() = eMesh(Vector{SVector{3,Float64}}(), Vector{SVector{3,Int64}}(), nothing, nothing)
    eMesh{Nothing,Tet}() = eMesh(Vector{SVector{3,Float64}}(), nothing, Vector{SVector{4,Int64}}(), Vector{Float64}())
    eMesh{Tri,Tet}() = eMesh(Vector{SVector{3,Float64}}(), Vector{SVector{3,Int64}}(), Vector{SVector{4,Int64}}(), Vector{Float64}())
end

as_tet_eMesh(e_mesh::eMesh{Tri,Tet}) = eMesh(e_mesh.point, nothing, e_mesh.tet, e_mesh.ϵ)
as_tri_eMesh(e_mesh::eMesh{Tri,Tet}) = eMesh(e_mesh.point, e_mesh.tri, nothing, nothing)
as_tri_eMesh(e_mesh::eMesh{Tri,Nothing}) = deepcopy(e_mesh)

vertex_pos_for_tri_ind(eM::eMesh{Tri,T2}, k::Int64) where {T2} = eM.point[eM.tri[k]]
vertex_pos_for_tet_ind(eM::eMesh{T1,Tet}, k::Int64) where {T1} = eM.point[eM.tet[k]]

n_point(eM::eMesh) = length(eM.point)
n_tri(eM::eMesh) = length(eM.tri)
n_tet(eM::eMesh) = length(eM.tet)
get_tri(eM::eMesh) = eM.tri
get_tet(eM::eMesh) = eM.tet
get_point(eM::eMesh) = eM.point

function Base.empty!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    empty!(e_mesh.point)
    (T1 == Nothing) || empty!(e_mesh.tri)
    (T2 == Nothing) || empty!(e_mesh.tet)
    (T2 == Nothing) || empty!(e_mesh.ϵ)
    return nothing
end

function Base.isempty(e_mesh::eMesh{T1,T2}) where {T1,T2}
    is_emp = isempty(e_mesh.point)
    if T1 == Tri
        is_emp = is_emp || isempty(e_mesh.tri)
    end
    if T2 == Tet
        is_emp = is_emp || isempty(e_mesh.tet)
        is_emp = is_emp || isempty(e_mesh.ϵ)
    end
    return is_emp
end

function Base.append!(eM_1::eMesh{T1,T2}, eM_2::eMesh{T1,T2}) where {T1,T2}
    n_1 = n_point(eM_1)
    if T1 != Nothing
        for k = eM_2.tri
            push!(eM_1.tri, k .+ n_1)
        end
    end
    if T2 != Nothing
        for k = eM_2.tet
            push!(eM_1.tet, k .+ n_1)
        end
        append!(eM_1.ϵ, eM_2.ϵ)
    end
    append!(eM_1.point, eM_2.point)
    return nothing
end

### VERIFICATION

function verify_mesh(eM::eMesh{T1,T2}) where {T1,T2}
    verify_mesh_triangles(eM)
    verify_mesh_tets(eM)
    return nothing
end

verify_mesh_tets(eM::eMesh{T1,Nothing}) where {T1} = nothing
function verify_mesh_tets(eM::eMesh{T1,Tet}) where {T1}
    length_tet = n_tet(eM)
    length_point = n_point(eM)
    length_ϵ = length(eM.ϵ)
    for k = 1:length_tet
        iΔ = eM.tet[k]
        all(1 .<= iΔ .<= length_point) || error("index in tet with sides $(iΔ) not within points")
    end
    (length_ϵ == length_point) || error("number of points not equal to number or ϵ")
    return nothing
end

verify_mesh_triangles(eM::eMesh{Nothing,T2}) where {T2} = nothing
function verify_mesh_triangles(eM::eMesh{Tri,T2}) where {T2}
    length_tri = n_tri(eM)
    length_point = n_point(eM)
    for k = 1:length_tri
        iΔ = eM.tri[k]
        all(1 .<= iΔ .<= length_point) || error("index in triangle with sides $(iΔ) not within points")
    end
    return nothing
end

### MESH MANIPULATION

function dh_transform_mesh!(e_mesh::eMesh{T1,T2}, dh::basic_dh{Float64}) where {T1,T2}
    point = e_mesh.point
    for k = 1:n_point(e_mesh)
        point[k] = dh_vector_mul(dh, point[k])
    end
    return nothing
end

function scale!(e_mesh::eMesh, r::Union{Float64,SVector{3,Float64}})  # TODO: add this functionality to basic_dh constructor
    r = ones(SVector{3,Float64}) .* r
    sv_33 = SMatrix{3,3,Float64,9}(r[1], 0.0, 0.0, 0.0, r[2], 0.0, 0.0, 0.0, r[3])
    dh_transform_mesh!(e_mesh, basic_dh(sv_33))
end

function crop_mesh(e_mesh::eMesh{Tri,Nothing}, n̂::SVector{3,Float64}, d::Float64, is_hard::Bool=false) # where {T2}
    Base.depwarn("This function is depricated, call crop_mesh(e_mesh, plane::SMatrix{1,4,Float64,4}) instead.", :crop_mesh)
    return crop_mesh(e_mesh, SMatrix{4,1,Float64,4}(n̂..., d), is_hard)
end

function crop_mesh(e_mesh::eMesh{Tri,Nothing}, plane::SMatrix{1,4,Float64,4}, is_hard::Bool=false) # where {T2}
    n̂ = unPad(plane)
    d = plane[4]
    e_mesh = deepcopy(e_mesh)
    mesh_repair!(e_mesh)
    if !isempty(e_mesh)
        p = get_point(e_mesh)
        i = get_tri(e_mesh)
        n_faces_orig = length(i)
        bool_delete = falses(n_faces_orig)
        area_ = [area(p[k]) for k = i]
        min_area = minimum(area_)
        min_area_tol = min_area / 100
        for k = 1:n_faces_orig
            recursivly_rotate!(i, p, bool_delete, k, n̂, d, 1, is_hard)
        end
        area_ = [area(p[k]) for k = i]
        i_small = findall(area_ .<= min_area_tol)
        i_delete = findall(bool_delete)
        append!(i_delete, i_small)
        i_delete = sort(unique(i_delete))
        deleteat!(i, i_delete)
        e_mesh_new = eMesh(p, i, nothing, nothing)
        mesh_repair!(e_mesh_new)
        return e_mesh_new
    end
    return eMesh{Tri,Nothing}()
end

function invert!(eM::eMesh{Tri,Nothing})
    for k = 1:n_tri(eM)
        tri_k = eM.tri[k]
        eM.tri[k] = SVector{3,Float64}(tri_k[3], tri_k[2], tri_k[1])
    end
end

### MESH REPAIR
rekey!(v::Nothing, i::Vector{Int64}) = nothing
rekey!(v::Vector{SVector{N,Int64}}, i::Vector{Int64}) where {N} = replace!(x -> i[x], v)

function mesh_inplace_rekey(e_mesh::eMesh{T1,T2}) where {T1,T2}
    e_mesh_copy = deepcopy(e_mesh)
    mesh_inplace_rekey!(e_mesh_copy)
    return e_mesh_copy
end

function mesh_inplace_rekey!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    # Expresses the mesh using the minumum number of points

    function shortest_side(e_mesh::eMesh)
        shortest_side(v::Nothing) = Inf
        shortest_side(ind::Vector{SVector{N,Int64}}) where {N} =  minimum([shortest_side(e_mesh.point[k]) for k = ind])
        function shortest_side(v::SVector{N,SVector{3,Float64}}) where {N}
            d_min = Inf
            for k = 1:N
                for kk = 1:(k-1)
                    d_min = min(d_min, norm(v[k] - v[kk]))
                end
            end
            return d_min
        end

        return min(shortest_side(e_mesh.tri), shortest_side(e_mesh.tet))
    end

    function inplace_rekey(point::Vector{SVector{3,Float64}}, min_side_length::Float64)
        balltree = BallTree(point; reorder = false)
        idxs = inrange(balltree, point, min_side_length * 0.499)  # Assume that all points closer to each other than hald a side length are duplicates
        return first.(sort!.(idxs))
    end

    if !isempty(e_mesh)  # reducing over an empty collection is not allowed
        min_side_length = shortest_side(e_mesh)
        new_key = inplace_rekey(e_mesh.point, min_side_length)
        rekey!(e_mesh.tri, new_key)
        rekey!(e_mesh.tet, new_key)
        mesh_remove_unused_points!(e_mesh)
        return new_key
    end
    return nothing
end

function mesh_remove_unused_points!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    truth_vector = falses(n_point(e_mesh))
    (T1 != Nothing) && (for k = e_mesh.tri; truth_vector[k] = true; end)
    (T2 != Nothing) && (for k = e_mesh.tet; truth_vector[k] = true; end)
    new_key = cumsum(truth_vector)
    rekey!(e_mesh.tri, new_key)
    rekey!(e_mesh.tet, new_key)
    ind_truth = findall(truth_vector)
    point_new = e_mesh.point[ind_truth]
    empty!(e_mesh.point)
    append!(e_mesh.point, point_new)
    if T2 != Nothing
        ϵ_new = e_mesh.ϵ[ind_truth]
        empty!(e_mesh.ϵ)
        append!(e_mesh.ϵ, ϵ_new)
    end
    return nothing
end

delete_triangles!(e_mesh::eMesh{Nothing,T2}) where {T2} = -9999
function delete_triangles!(e_mesh::eMesh{Tri,T2}) where {T2}
    ### make the first index lowest
    for k = 1:n_tri(e_mesh)
        iΔ = e_mesh.tri[k]
        min_index = minimum(iΔ)
        (min_index == iΔ[1]) || (iΔ = SVector(iΔ[3], iΔ[1], iΔ[2]))
        (min_index == iΔ[1]) || (iΔ = SVector(iΔ[3], iΔ[1], iΔ[2]))
    end

    # println("e_mesh.tri: ", e_mesh.tri)

    ### create a dictionary of key repition counts
    key_type = Vector{Tuple{Int64,SVector{3,Int64}}}
    dict_ind = Dict{SVector{3,Int64},key_type}()
    for k = 1:n_tri(e_mesh)
        iΔ = e_mesh.tri[k]
        sort_iΔ = sort(iΔ)
        !haskey(dict_ind, sort_iΔ) && (dict_ind[sort_iΔ] = key_type() )
        push!(dict_ind[sort_iΔ], (k, iΔ))
    end

    # println("dict_ind: ", dict_ind)

    ### delete duplicates
    i_delete = Vector{Int64}()
    for (key_k, val_k) = dict_ind
        # println("val_k: ", val_k)
        length_val_k = length(val_k)
        if length_val_k == 2
            val_k1 = val_k[1]
            val_k2 = val_k[2]
            (val_k1[2] == val_k2[2]) && error("non-opposing triangles")
            push!(i_delete, val_k1[1], val_k2[1])
        elseif 3 <= length_val_k
            # println("length_key_k: ", length_val_k)
            # println("val_k: ", val_k)
            error("something is wrong")
        end
    end
end
# function delete_triangles!(e_mesh::eMesh{Tri,T2}) where {T2}
#     sort!(e_mesh.tri, by = x -> sort(x))
#     n_tri_delete = 0
#     for k = n_tri(e_mesh):-1:2
#         if k <= n_tri(e_mesh)
#             vert_2 = e_mesh.tri[k]
#             vert_1 = e_mesh.tri[k - 1]
#             if sort(vert_1) == sort(vert_2)  # share same 3 indices
#                 n̂_2 = triangleNormal(e_mesh.point[vert_2])
#                 n̂_1 = triangleNormal(e_mesh.point[vert_1])
#                 if n̂_2 ≈ n̂_1  # triangle is a duplicate (normals point in same direction)
#                     deleteat!(e_mesh.tri, k - 1)  # delete this triangle
#                     n_tri_delete += 1
#                 elseif n̂_2 ≈ -n̂_1  # triangles oppose each other (normals point in different directions)
#                     is_check_this_k = false
#                     deleteat!(e_mesh.tri, k - 1)  # delete this triangle
#                     deleteat!(e_mesh.tri, k - 1)  # delete this triangle
#                     n_tri_delete += 2
#                 else
#                     println("n̂_1: ", n̂_1)
#                     println("n̂_2: ", n̂_2)
#                     error("something is wrong")
#                 end
#             end
#         end
#     end
#     return n_tri_delete
# end

function sub_div_mesh(eM_ico::eMesh{Tri,T2}, n_div::Int64) where {T2}
    n_end(n::Int64) = div((n + 1) * n, 2)
    n_start(n::Int64) = 1 + n_end(n - 1)

    function sub_div_triangle(p::SVector{3,SVector{3,Float64}}, n_div::Int64)
        function sub_div_triangle_vert_index(n_div::Int64)
            i_tri = Vector{SVector{3,Int64}}()
            for k = 1:n_div
                for kk = 0:(k - 1)
                    i1 = n_start(k) + kk
                    i2 = i1 + k
                    i3 = i2 + 1
                    push!(i_tri, SVector{3,Int64}(i1, i2, i3))
                end
                for kk = 0:(k - 2)
                    i1 = n_start(k) + kk
                    i2 = i1 + k + 1
                    i3 = i2 - k
                    push!(i_tri, SVector{3,Int64}(i1, i2, i3))
                end
            end
            return i_tri
        end

        function get_new_point(n_vert::Int64, n_div::Int64, p::SVector{3,SVector{3,Float64}})
            function find_layer_1(i_point::Int64)
                i_end_layer = 1
                i_layer = 1
                while i_end_layer < i_point
                    i_layer += 1
                    i_end_layer += i_layer
                end
                norm_extent_1 = ifelse(i_point == 1, 0.0, (i_end_layer - i_point) / (i_layer - 1) )
                return i_layer, norm_extent_1
            end

            i_layer, norm_extent_1 = find_layer_1(n_vert)
            ϕ_1 = (n_div - i_layer + 1) / n_div
            ϕ_2 = (1 - ϕ_1) * norm_extent_1
            ϕ_3 = 1 - ϕ_1 - ϕ_2
            ϕ = SVector{3,Float64}(ϕ_1, ϕ_2, ϕ_3)
            return sum(p .* ϕ)
        end

        point = Vector{SVector{3,Float64}}()
        i_tri_div = sub_div_triangle_vert_index(n_div)
        for k = 1:n_end(n_div + 1)
            push!(point, get_new_point(k, n_div, p))
        end
        return eMesh(point, i_tri_div, nothing, nothing)
    end

    eM_ico_div = eMesh{Tri,Nothing}()
    for k = 1:n_tri(eM_ico)
        p = vertex_pos_for_tri_ind(eM_ico, k)
        append!(eM_ico_div, sub_div_triangle(p, n_div))
    end
    mesh_repair!(eM_ico_div)
    return eM_ico_div
end

function mesh_repair!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    mesh_remove_unused_points!(e_mesh)
    mesh_inplace_rekey!(e_mesh)
    mesh_remove_unused_points!(e_mesh)
    return delete_triangles!(e_mesh)
end

### BASIC SHAPES
function output_eMesh_half_plane(plane_w::Float64=1.0, is_include_vis_sides::Bool=false)
    v1, v2, v3 = [SVector{3,Float64}(cos(theta), sin(theta), 0.0) for theta = (0.0, 2*pi/3, 4*pi/3)]
    v4 = SVector{3,Float64}(0.0, 0.0, -1.0) * plane_w
    point = [v1, v2, v3, v4]
    if is_include_vis_sides
        tri = [SVector{3,Int64}(1,2,3), SVector{3,Int64}(1,3,4),  SVector{3,Int64}(1,4,2), SVector{3,Int64}(2,4,3)]
    else
        tri = [SVector{3,Int64}(1,2,3)]
    end
    tet = [SVector{4,Int64}(4,1,2,3)]
    ϵ = [0.0, 0.0, 0.0, -plane_w]
    return eMesh(point, tri, tet, ϵ)
end

function output_eMesh_sphere(rad::Float64=1.0, n_div::Int64=4)
    function make_icosahedron()
        φ = Base.MathConstants.golden

        v = Vector{SVector{3,Float64}}()
        for s_1 = (-1.0, +1.0)
            for s_2 = (-1.0, +1.0)
                push!(v, SVector{3,Float64}(    0.0,     s_1, φ * s_2) )
                push!(v, SVector{3,Float64}(    s_1, φ * s_2,     0.0) )
                push!(v, SVector{3,Float64}(φ * s_2,     0.0,     s_1) )
            end
        end

        n_ico_vert = 12
        d = zeros(n_ico_vert, n_ico_vert)
        for k = 1:n_ico_vert
            for kk = 1:n_ico_vert
                d[k, kk] = norm(v[k] - v[kk])
            end
        end

        v_face_vert = Vector{SVector{3,Int64}}()
        b = d .== 2.0
        for i1 = 1:n_ico_vert
            for i2 = 1:n_ico_vert
                for i3 = 1:n_ico_vert
                    if (i1 < i2 < i3) && b[i1, i2] && b[i2, i3] && b[i1, i3]
                        p1 = v[i1]
                        p2 = v[i2]
                        p3 = v[i3]
                        n = triangleNormal(p1, p2, p3)
                        c = normalize(centroid(p1, p2, p3))
                        i_face = ifelse(n ≈ c, SVector{3,Int64}(i1, i2, i3), SVector{3,Int64}(i1, i3, i2))
                        push!(v_face_vert, i_face)
                    end
                end
            end
        end

        push!(v, SVector{3,Float64}(0,0,0))
        v_tet = Vector{SVector{4,Int64}}()
        for k = 1:length(v_face_vert)
            push!(v_tet, SVector{4,Int64}(13, v_face_vert[k]...))
        end
        ϵ = vcat(zeros(n_ico_vert), -1.0)
        return eMesh(v, v_face_vert, v_tet, ϵ)
    end

    function volumize_about(eM::eMesh{Tri,Nothing})
        i_tet = Vector{SVector{4,Int64}}()
        n_vert = n_point(eM)
        n_center = n_vert + 1
        for k = 1:n_tri(eM)
            push!(i_tet, SVector{4,Int64}(n_center, eM.tri[k]...))
        end
        ϵ = vcat(zeros(n_vert), -1.0)
        point = deepcopy(eM.point)
        push!(point, zeros(SVector{3,Float64}))
        return eMesh(point, eM.tri, i_tet, ϵ)
    end

    function project_to_sphere!(eM::eMesh, rad::Float64=1.0)
        for k = 1:n_point(eM)
            eM.point[k] = normalize(eM.point[k]) * rad
        end
        return nothing
    end

    eM_ico = make_icosahedron()
    eM_ico_div = sub_div_mesh(eM_ico, n_div)
    project_to_sphere!(eM_ico_div, rad)
    return volumize_about(eM_ico_div)
end

function make_cone_points(rⁱ::Float64, rᵒ::Float64, h::Float64, k_slice::Int64, tot_slice::Int64)
    polar_point(r::Float64, θ::Float64, z::Float64) = SVector{3,Float64}(r * cos(θ), r * sin(θ), z)

    @assert(0 < rⁱ < rᵒ)
    @assert(0 < h)
    @assert(1 <= k_slice <= tot_slice)

    v = Vector{SVector{3,Float64}}()
    θ_space = LinRange{Float64}(0.0, 2*pi, tot_slice + 1)
    θ_1 = θ_space[k_slice]
    θ_2 = θ_space[k_slice + 1]
    for k = 1:8
        θ = ifelse(mod1(k, 4) <= 2, θ_1, θ_2)
        r = ifelse(isodd(k), rⁱ, rᵒ)
        z = ifelse(k <= 4, -0.5, +0.5) * h
        push!(v, polar_point(r, θ, z))
    end
    push!(v, polar_point((rⁱ + rᵒ) / 2, (θ_1 + θ_2) / 2, 0.0))
end

function output_box_ind()
    oriented_box_faces = (
        SVector{4, Int64}(1,3,5,7),  # -x
        SVector{4, Int64}(2,6,4,8),  # +x
        SVector{4, Int64}(1,5,2,6),  # -y
        SVector{4, Int64}(3,4,7,8),  # +y
        SVector{4, Int64}(1,2,3,4),  # -z
        SVector{4, Int64}(5,7,6,8),  # +z
    )
    ϵ = zeros(Float64,8)
    push!(ϵ, -1.0)
    tri = Vector{SVector{3,Int64}}()
    tet = Vector{SVector{4,Int64}}()
    for k = 1:6
        bf_k = oriented_box_faces[k]
        push!(tri, bf_k[SVector{3,Int64}(1,3,4)])
        push!(tri, bf_k[SVector{3,Int64}(1,4,2)])
    end
    for k = 1:length(tri)
        tri_k = tri[k]
        push!(tet, SVector{4,Int64}(9, tri_k[1], tri_k[2], tri_k[3]))
    end
    return tri, tet, ϵ
end

function output_eMesh_box(r::Union{Float64,SVector{3,Float64}}=1.0, c::SVector{3,Float64}=zeros(SVector{3,Float64}))
    point = [
        SVector{3,Float64}(-1,-1,-1),
        SVector{3,Float64}(+1,-1,-1),
        SVector{3,Float64}(-1,+1,-1),
        SVector{3,Float64}(+1,+1,-1),
        SVector{3,Float64}(-1,-1,+1),
        SVector{3,Float64}(+1,-1,+1),
        SVector{3,Float64}(-1,+1,+1),
        SVector{3,Float64}(+1,+1,+1),
        SVector{3,Float64}( 0, 0, 0),
    ]
    tri, tet, ϵ = output_box_ind()
    e_mesh = eMesh(point, tri, tet, ϵ)
    scale!(e_mesh, r)
    dh_transform_mesh!(e_mesh, basic_dh(c))
    return e_mesh
end

function output_eMesh_slice(rⁱ::Float64, rᵒ::Float64, h::Float64, k_slice::Int64, tot_slice::Int64)
    point = make_cone_points(rⁱ, rᵒ, h, k_slice, tot_slice)
    tri, tet, ϵ = output_box_ind()
    return eMesh(point, tri, tet, ϵ)
end

function output_eMesh_hole(rⁱ::Float64, rᵒ::Float64, h::Float64, tot_slice::Int64)
    eM_cone = eMesh{Tri,Tet}()
    for k = 1:tot_slice
        append!(eM_cone, output_eMesh_slice(rⁱ, rᵒ, h, k, tot_slice))
    end
    return eM_cone
end

# function eMesh{Tri,Nothing}()
#     point = Vector{SVector{3,Float64}}()
#     tri = Vector{SVector{3,Int64}}()
#     return eMesh(point, tri, nothing, nothing)
# end

# function eMesh{Nothing,Tet}()
#     point = Vector{SVector{3,Float64}}()
#     tri = nothing
#     tet = Vector{SVector{4,Int64}}()
#     ϵ = Vector{Float64}()
#     return eMesh(point, nothing, tet, ϵ)
# end

# function eMesh{Tri,Tet}()
#     point = Vector{SVector{3,Float64}}()
#     tri = Vector{SVector{3,Int64}}()
#     tet = Vector{SVector{4,Int64}}()
#     ϵ = Vector{Float64}()
#     return eMesh(point, tri, tet, ϵ)
# end
