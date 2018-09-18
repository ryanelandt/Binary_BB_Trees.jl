is_good_orientation(v) = (v[1] == 1) || (v[1] == 8)
rotateIndexAboutZ(v) = v[[3, 1, 4, 2, 7, 5, 8, 6]]
cycleArray!(v) = pushfirst!(v, pop!(v))
function makeMinFirst!(ix)
    cycleArray!(ix)
    (findmin(ix)[2] == 1) && (return nothing)
    makeMinFirst!(ix)
end
function outputOrientedBoxFaces()
    return [SVector{4, Int64}(1,3,5,7),  # -x
            SVector{4, Int64}(2,6,4,8),  # +x
            SVector{4, Int64}(1,5,2,6),  # -y
            SVector{4, Int64}(3,4,7,8),  # +y
            SVector{4, Int64}(1,2,3,4),  # -z
            SVector{4, Int64}(5,7,6,8)]  # +z
end
function outputDividedCube()
    a6 = Vector{SVector{4,Int64}}()
    i_face_all = outputOrientedBoxFaces()
    for k_cube_side = 1:6
        i_face_k = Vector{Int64}(i_face_all[k_cube_side])
        v = vcat(i_face_k, i_face_k .+ 8)
        !is_good_orientation(v) && (v = rotateIndexAboutZ(v))
        !is_good_orientation(v) && (v = rotateIndexAboutZ(v))
        !is_good_orientation(v) && (v = rotateIndexAboutZ(v))
        !is_good_orientation(v) && error("should be centered now")
        magic_num = v[1]
        for s_xyz = ([2, 4, 8, 6], [4, 3, 7, 8], [5, 6, 8, 7])
            ix = v[s_xyz]
            makeMinFirst!(ix)
            push!(a6, SVector{4,Int64}(magic_num, ix[1], ix[2], ix[3]))
            push!(a6, SVector{4,Int64}(magic_num, ix[1], ix[3], ix[4]))
        end
    end
    return a6
end
function determineUniquePoints(points_16)
    i_unique = 0
    dict_points = Dict{SVector{3,Float64},Int64}()
    for (k, point_k) = enumerate(points_16)
        if !haskey(dict_points, point_k)
            i_unique += 1
            dict_points[point_k] = i_unique
        end
    end
    return dict_points
end
function outputBoxTriInd()
    v_tri = Vector{SVector{3,Int64}}()
    [twoTriangles!(v_tri, k) for k = outputOrientedBoxFaces()]
    return v_tri
end
box_indices() = [SVector(Tuple(k_cart)) for k_cart = CartesianIndices((2,2,2))[:]]
basicBoxPoints() = [SVector{3,Float64}(ifelse.(Tuple(k).==1,-1.0,+1.0)) for k = box_indices()]
function sizeCenterBoxPoints(scale::SVector{3,Float64}=SVector{3,Float64}(1,1,1), center::SVector{3,Float64}=SVector{3,Float64}(0,0,0))
    v = basicBoxPoints()
    for k = eachindex(v)
        v[k] = v[k] .* scale .+ center
    end
    return v
end
function twoTriangles!(v_tri, i_in)
    push!(v_tri, i_in[SVector{3, Int64}(1,3,4)])
    push!(v_tri, i_in[SVector{3, Int64}(1,4,2)])
    return nothing
end
function outputBoxVolMesh(rad_box, foam_depth, center)
    all(0.01 .<= foam_depth .<= 1.0) || error("foam depth needs to be between 0.01 and 1.0 of box radius")
    rad_inner = rad_box .* (1 .- foam_depth)
    points_16 = vcat(sizeCenterBoxPoints(rad_box, center), sizeCenterBoxPoints(rad_inner, center))
    dict_points = determineUniquePoints(points_16)
    vec_nondegenerate_tet = Vector{SVector{4,Int64}}()
    for ind_k = outputDividedCube()
        ss_k = SVector{4,Int64}([dict_points[k] for k = points_16[ind_k]])
        (length(unique(ss_k)) == 4) && push!(vec_nondegenerate_tet, ss_k)
    end
    points_12 = Vector{SVector{3,Float64}}(undef, length(dict_points))
    for (key_k, k) = dict_points
        points_12[k] = key_k
    end
    return points_12, outputBoxTriInd(), vec_nondegenerate_tet
end


# # @inline basicBoxPoints() = [SVector{3,Float64}(ifelse.(Tuple(k).==1,-1.0,+1.0)) for k = CartesianIndices((2,2,2))[:]]
# box_indices() = [SVector(Tuple(k_cart)) for k_cart = CartesianIndices((2,2,2))[:]]
# basicBoxPoints() = [SVector{3,Float64}(ifelse.(Tuple(k).==1,-1.0,+1.0)) for k = box_indices()]
# function cubeToFiveTets(i_in)
#     i1 = SVector{4, Int64}ze about equilibrium points; in these cases u0 is a particular one that maintains the equilibrium. For non-equilibrium linearizations it can be anything. Linearization about a non-equilibrium state may occur during non-linear optimization (e.g. as part of a gradient descent scheme).(1,2,4,6)
#     i2 = SVector{4, Int64}(1,4,3,7)
#     i3 = SVector{4, Int64}(1,6,7,5)
#     i4 = SVector{4, Int64}(4,7,6,8)
#     i5 = SVector{4, Int64}(1,4,7,6)  # bigger one
#     return i_in[i1], i_in[i2], i_in[i3], i_in[i4], i_in[i5]
# end
# function twoTriangles!(v_tri, i_in)
#     i1 = SVector{3, Int64}(1,3,4)
#     i2 = SVector{3, Int64}(1,4,2)
#     push!(v_tri, i_in[i1])
#     push!(v_tri, i_in[i2])
#     return nothing
# end
# function outputOrientedBoxFaces()
#     return [SVector{4, Int64}(1,3,5,7),  # -x
#             SVector{4, Int64}(2,6,4,8),  # +x
#             SVector{4, Int64}(1,5,2,6),  # -y
#             SVector{4, Int64}(3,4,7,8),  # +y
#             SVector{4, Int64}(1,2,3,4),  # -z
#             SVector{4, Int64}(5,7,6,8)]  # +z
# end
# function sizeCenterBoxPoints(scale::SVector{3,Float64}=SVector{3,Float64}(1,1,1), center::SVector{3,Float64}=SVector{3,Float64}(0,0,0))
#     v = basicBoxPoints()
#     for k = eachindex(v)
#         v[k] = v[k] .* scale .+ center
#     end
#     return v
# end
# function newTruth(is_reduce::SVector{3,Bool})
#     i_9_16 = collect(reshape(9:16, 2, 2, 2))
#     for (k, the_ind) = enumerate(CartesianIndices(i_9_16))
#         i1, i2, i3 = Tuple(the_ind)
#         for j = 1:3
#             mm = MVector{3,Int64}(Tuple(the_ind))
#             mm[j] = 3 - mm[j]
#             is_reduce[j] && (i_9_16[k] = min(i_9_16[k], i_9_16[mm[1], mm[2], mm[3]]))
#         end
#     end
#     return vcat(collect(1:8), i_9_16[:])
# end
# function outputBoxTet30()
#     vec_mesh_tet = Vector{SVector{4,Int64}}()
#     i_face = outputOrientedBoxFaces()
#     for (k, i_face_k) = enumerate(i_face)
#         vec_mesh_tet_new = [cubeToFiveTets(vcat(i_face_k, i_face_k .+ 8))...]
#         append!(vec_mesh_tet, vec_mesh_tet_new)
#     end
#     return vec_mesh_tet
# end
# function outputBoxTriInd()
#     v_tri = Vector{SVector{3,Int64}}()
#     [twoTriangles!(v_tri, k) for k = outputOrientedBoxFaces()]
#     return v_tri
# end
# function boxMesh(rad_box::SVector{3,Float64}=SVector{3,Float64}(1,1,1), foam_depth::SVector{3,Float64}=SVector{3,Float64}(0,0,0), center::SVector{3,Float64}=SVector{3,Float64}(0,0,0))
#     rad_inner = rad_box .* (1.0 .- foam_depth)
#     points_30 = vcat(sizeCenterBoxPoints(rad_box, center), sizeCenterBoxPoints(rad_inner, center))
#     is_reduce = foam_depth .== 1.0
#     # find mapings between reduced and master cube
#     i_vert_reduce = newTruth(is_reduce)
#     unique_vert_reduce = sort(unique(i_vert_reduce))
#     ind_orig_2_ind_new = [findfirst(unique_vert_reduce .== k) for k = i_vert_reduce]
#     # purge degenerate tets
#     vec_ind_tet = [ind_orig_2_ind_new[k] for k = outputBoxTet30()]
#     vec_ind_tet = vec_ind_tet[[length(unique(k))==4 for k = vec_ind_tet]]
#     # Make Mesh
#     return SurfaceVolumeMesh(points_30[unique_vert_reduce], outputBoxTriInd(), vec_ind_tet)
# end
