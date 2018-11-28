
abstract type Tri end
abstract type Tet end

struct eMesh{T1<:Union{Nothing,Tri},T2<:Union{Nothing,Tet}}
    point::Vector{SVector{3,Float64}}
    tri::Union{Nothing,Vector{SVector{3,Int64}}}
    tet::Union{Nothing,Vector{SVector{4,Int64}}}
    function eMesh(point::Vector{SVector{3,Float64}},
                          tri::Union{Nothing,Vector{SVector{3,Int64}}},
                          tet::Union{Nothing,Vector{SVector{4,Int64}}}=nothing)
        #
        T1_ = ifelse(tri == nothing, Nothing, Tri)
        T2_ = ifelse(tet == nothing, Nothing, Tet)
        T1_ == T2_ == Nothing && error("a whole lot of nothing")
        return new{T1_,T2_}(point, tri, tet)
    end
    function eMesh(hm::HomogenousMesh, tet::Union{Nothing,Vector{SVector{4,Int64}}}=nothing)
        point = get_h_mesh_vertices(hm)
        tri = get_h_mesh_faces(hm)
        return eMesh(point, tri, tet)
    end
end

n_points(eM::eMesh) = length(eM.point)
n_tri(eM::eMesh) = length(eM.tri)
n_tet(eM::eMesh) = length(eM.tet)

function Base.empty!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    (T1 == Nothing) || empty!(e_mesh.tri)
    (T2 == Nothing) || empty!(e_mesh.tet)
    empty!(e_mesh.point)
    return nothing
end

function Base.append!(eM_1::eMesh{T1,T2}, eM_2::eMesh{T1,T2}) where {T1,T2}
    n_1 = n_points(eM_1)
    if T1 != Nothing
        for k = eM_2.tri
            push!(eM_1.tri, k .+ n_1)
        end
    end
    if T2 != Nothing
        for k = eM_2.tet
            push!(eM_1.tet, k .+ n_1)
        end
    end
    append!(eM_1.point, eM_2.point)
    return nothing
end

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
    truth_vector = falses(n_points(e_mesh))
    (T1 != Nothing) && (for k = eM_2.tri; truth_vector[k] = true; end)
    (T2 != Nothing) && (for k = eM_2.tet; truth_vector[k] = true; end)
    new_key = cumsum(truth_vector)
    rekey!(e_mesh.tri, new_key)
    rekey!(e_mesh.tet, new_key)
    point_new = e_mesh.point[findall(truth_vector)]
    empty!(e_mesh.point)
    append!(e_mesh.point, point_new)
    return nothing
end
