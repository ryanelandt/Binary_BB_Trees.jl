
abstract type Tri end
abstract type Tet end

struct eMesh{T1<:Union{Nothing,Tri},T2<:Union{Nothing,Tet}}
    point::Vector{SVector{3,Float64}}
    tri::Union{Nothing,Vector{SVector{3,Int64}}}
    tet::Union{Nothing,Vector{SVector{4,Int64}}}
    ϵ::Union{Nothing,Vector{Float64}}
    function eMesh( point::Vector{SVector{3,Float64}},
                    tri::Union{Nothing,Vector{SVector{3,Int64}}},
                    tet::Union{Nothing,Vector{SVector{4,Int64}}},
                    ϵ::Union{Nothing,Vector{Float64}})

        T1_ = ifelse(tri == nothing, Nothing, Tri)
        T2_ = ifelse(tet == nothing, Nothing, Tet)
        if T2_ == Tet
            @assert(isa(ϵ, Vector{Float64}))
            @assert(length(ϵ) == length(point), "length(ϵ) = $(length(ϵ)) but length(point) = $(length(point))")
            (length(ϵ) != 0) && @assert(maximum(ϵ) == 0.0, "strain must be zero on the surface of the volume mesh")
            for k = 1:length(tet)
                (0.0 < volume(point[tet[k]])) || error("something is wrong")
            end
        end
        T1_ == T2_ == Nothing && error("a whole lot of nothing")
        return new{T1_,T2_}(point, tri, tet, ϵ)
    end
    function eMesh(hm::HomogenousMesh, tet::Union{Nothing,Vector{SVector{4,Int64}}},
            ϵ::Union{Nothing,Vector{Float64}})

        point = get_h_mesh_vertices(hm)
        tri = get_h_mesh_faces(hm)
        return eMesh(point, tri, tet, ϵ)
    end
    function eMesh{Tri,Nothing}()
        point = Vector{SVector{3,Float64}}()
        tri = Vector{SVector{3,Int64}}()
        return eMesh(point, tri, nothing, nothing)
    end
    function eMesh{Nothing,Tet}()
        point = Vector{SVector{3,Float64}}()
        tri = nothing
        tet = Vector{SVector{4,Int64}}()
        ϵ = Vector{Float64}()
        return eMesh(point, nothing, tet, ϵ)
    end
    function eMesh{Tri,Tet}()
        point = Vector{SVector{3,Float64}}()
        tri = Vector{SVector{3,Int64}}()
        tet = Vector{SVector{4,Int64}}()
        ϵ = Vector{Float64}()
        return eMesh(point, tri, tet, ϵ)
    end
end

as_tet_eMesh(e_mesh::eMesh{Tri,Tet}) = eMesh(e_mesh.point, nothing, e_mesh.tet, e_mesh.ϵ)
as_tri_eMesh(e_mesh::eMesh{Tri,Tet}) = eMesh(e_mesh.point, e_mesh.tri, nothing, nothing)

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
    empty!(e_mesh.ϵ)
    return nothing
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

### MESH MANIPULATION

function dh_transform_mesh!(e_mesh::eMesh{T1,T2}, dh::basic_dh{Float64}) where {T1,T2}
    point = e_mesh.point
    for k = 1:n_point(e_mesh)
        point[k] = dh_vector_mul(dh, point[k])
    end
    return nothing
end

function scale!(e_mesh::eMesh, r::Union{Float64,SVector{3,Float64}})
    r = ones(SVector{3,Float64}) .* r
    sv_33 = SMatrix{3,3,Float64,9}(r[1], 0.0, 0.0, 0.0, r[2], 0.0, 0.0, 0.0, r[3])
    dh_transform_mesh!(e_mesh, basic_dh(sv_33))
end

function crop_mesh(e_mesh::eMesh{Tri,T2}, n̂::SVector{3,Float64}, d::Float64, is_hard::Bool=false) where {T2}
    (T2 == Nothing) || @warn "tetrahedron clipping NOT handled correctly"
    m = deepcopy(e_mesh)
    p = get_point(m)
    i = get_tri(m)
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
    ϵ = zeros(Float64, length(p))

    if T2 == Nothing
        e_mesh_new = eMesh(p, i, nothing, nothing)
    else
        e_mesh_new = eMesh(p, i, e_mesh.tet, ϵ)
    end

    mesh_repair!(e_mesh_new)
    return e_mesh_new
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

    min_side_length = shortest_side(e_mesh)
    new_key = inplace_rekey(e_mesh.point, min_side_length)
    rekey!(e_mesh.tri, new_key)
    rekey!(e_mesh.tet, new_key)
    return nothing
end

function mesh_remove_unused_points!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    truth_vector = falses(n_point(e_mesh))
    (T1 != Nothing) && (for k = e_mesh.tri; truth_vector[k] = true; end)
    (T2 != Nothing) && (for k = e_mesh.tet; truth_vector[k] = true; end)
    new_key = cumsum(truth_vector)
    rekey!(e_mesh.tri, new_key)
    rekey!(e_mesh.tet, new_key)

    point_new = e_mesh.point[findall(truth_vector)]
    empty!(e_mesh.point)
    append!(e_mesh.point, point_new)

    if T2 != Nothing
        ϵ_new = e_mesh.ϵ[findall(truth_vector)]
        empty!(e_mesh.ϵ)
        append!(e_mesh.ϵ, ϵ_new)
    end
    return nothing
end

delete_triangles!(e_mesh::eMesh{Nothing,T2}) where {T2} = -9999
function delete_triangles!(e_mesh::eMesh{Tri,T2}) where {T2}
    sort!(e_mesh.tri, by = x -> sort(x))
    n_tri_delete = 0
    for k = n_tri(e_mesh):-1:2
        if k <= n_tri(e_mesh)
            vert_2 = e_mesh.tri[k]
            vert_1 = e_mesh.tri[k - 1]
            if sort(vert_1) == sort(vert_2)  # share same 3 indices
                n̂_2 = triangleNormal(e_mesh.point[vert_2])
                n̂_1 = triangleNormal(e_mesh.point[vert_1])
                if n̂_2 ≈ n̂_1  # triangle is a duplicate (normals point in same direction)
                    deleteat!(e_mesh.tri, k - 1)  # delete this triangle
                    n_tri_delete += 1
                elseif n̂_2 ≈ -n̂_1  # triangles oppose each other (normals point in different directions)
                    is_check_this_k = false
                    deleteat!(e_mesh.tri, k - 1)  # delete this triangle
                    deleteat!(e_mesh.tri, k - 1)  # delete this triangle
                    n_tri_delete += 2
                else
                    println("n̂_1: ", n̂_1)
                    println("n̂_2: ", n̂_2)
                    error("something is wrong")
                end
            end
        end
    end
    return n_tri_delete
end

function mesh_repair!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    mesh_inplace_rekey!(e_mesh)
    mesh_remove_unused_points!(e_mesh)
    return delete_triangles!(e_mesh)
end

### BASIC SHAPES
function output_eMesh_half_plane(plane_w::Float64=1.0)
    v1, v2, v3 = [SVector{3,Float64}(cos(theta), sin(theta), 0.0) for theta = (0.0, 2*pi/3, 4*pi/3)]
    v4 = SVector{3,Float64}(0.0, 0.0, -1.0) * plane_w
    point = [v1, v2, v3, v4]
    tri = [SVector{3,Int64}(1,2,3)]
    tet = [SVector{4,Int64}(4,1,2,3)]
    ϵ = [0.0, 0.0, 0.0, -plane_w]
    return eMesh(point, tri, tet, ϵ)
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
