const SV3 = SVector{3,Float64}  # TODO: remove this

@inline basicBoxPoints() = [SV3(ifelse.(Tuple(k).==1,-1.0,+1.0)) for k = CartesianIndices((2,2,2))[:]]
function cubeToFiveTets(i_in)
    i1 = SVector{4, Int64}(1,2,4,6)
    i2 = SVector{4, Int64}(1,4,3,7)
    i3 = SVector{4, Int64}(1,6,7,5)
    i4 = SVector{4, Int64}(4,7,6,8)
    i5 = SVector{4, Int64}(1,4,7,6)
    return i_in[i1], i_in[i2], i_in[i3], i_in[i4], i_in[i5]
end
function twoTriangles!(v_tri, i_in)
    i1 = SVector{3, Int64}(1,3,4)
    i2 = SVector{3, Int64}(1,4,3)
    push!(v_tri, i_in[i1])
    push!(v_tri, i_in[i2])
    return nothing
end
function outputOrientedBoxFaces()
    return [SVector{4, Int64}(1,3,5,7),  # -x
            SVector{4, Int64}(2,6,4,8),  # +x
            SVector{4, Int64}(1,5,2,6),  # -y
            SVector{4, Int64}(3,4,7,8),  # +y
            SVector{4, Int64}(1,2,3,4),  # -z
            SVector{4, Int64}(5,7,6,8)]  # +z
end
function sizeCenterBoxPoints(scale::SVector{3,Float64}=SV3(1,1,1), center::SVector{3,Float64}=SV3(0,0,0))
    v = basicBoxPoints()
    for k = eachindex(v)
        v[k] = v[k] .* scale .+ center
    end
    return v
end
function newTruth(is_reduce::SVector{3,Bool})
    i_9_16 = collect(reshape(9:16, 2, 2, 2))
    for (k, the_ind) = enumerate(CartesianIndices(i_9_16))
        i1, i2, i3 = Tuple(the_ind)
        for j = 1:3
            mm = MVector{3,Int64}(Tuple(the_ind))
            mm[j] = 3 - mm[j]
            is_reduce[j] && (i_9_16[k] = min(i_9_16[k], i_9_16[mm[1], mm[2], mm[3]]))
        end
    end
    return vcat(collect(1:8), i_9_16[:])
end
function outputBoxTet30()
    vec_mesh_tet = Vector{SVector{4,Int64}}()
    i_face = outputOrientedBoxFaces()
    for (k, i_face_k) = enumerate(i_face)
        vec_mesh_tet_new = [cubeToFiveTets(vcat(i_face_k, i_face_k .+ 8))...]
        append!(vec_mesh_tet, vec_mesh_tet_new)
    end
    return vec_mesh_tet
end
function outputBoxTriInd()
    v_tri = Vector{SVector{3,Int64}}()
    [twoTriangles!(v_tri, k) for k = outputOrientedBoxFaces()]
    return v_tri
end
function boxMesh(rad_box::SV3=SV3(1,1,1), foam_depth::SV3=SV3(0,0,0), center::SV3=SV3(0,0,0))
    rad_inner = rad_box .* (1.0 .- foam_depth)
    points_30 = vcat(sizeCenterBoxPoints(rad_box, center), sizeCenterBoxPoints(rad_inner, center))
    is_reduce = foam_depth .== 1.0
    # find mapings between reduced and master cube
    i_vert_reduce = newTruth(is_reduce)
    unique_vert_reduce = sort(unique(i_vert_reduce))
    ind_orig_2_ind_new = [findfirst(unique_vert_reduce .== k) for k = i_vert_reduce]
    # purge degenerate tets
    vec_ind_tet = [ind_orig_2_ind_new[k] for k = outputBoxTet30()]
    vec_ind_tet = vec_ind_tet[[length(unique(k))==4 for k = vec_ind_tet]]
    # Make Mesh
    return SurfaceVolumeMesh(points_30[unique_vert_reduce], outputBoxTriInd(), vec_ind_tet)
end
