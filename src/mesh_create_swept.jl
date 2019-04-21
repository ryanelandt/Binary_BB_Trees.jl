
function f_swept_circle(r::Float64, θ::Float64)
    n̂_1 = SVector{3,Float64}( cos(θ), sin(θ), 0.0)
    n̂_2 = SVector{3,Float64}(-sin(θ), cos(θ), 0.0)
    return r * n̂_1, n̂_1, n̂_2
end

function f_swept_triv(θ::Float64)
	n̂_1 = SVector(0.0, 0.0, -1.0)
	n̂_2 = SVector(0.0, 1.0, 0.0)
	return n̂_2 * θ, n̂_1, n̂_2
end

function add_rot_sym_segment!(eM, fun_gen, θ, ϕ, rad, is_open::NTuple{2,Bool})
    # x̂ is in the "radial" direction
    # ŷ is along the path

    p1, x̂1, ŷ1 = fun_gen(θ[1])
    p2, x̂2, ŷ2 = fun_gen(θ[2])
    p3 = (p1 + p2) * 0.5
    p4 = p1 + AngleAxis(ϕ[1], ŷ1...) * x̂1 * rad[1]
    p6 = p1 + AngleAxis(ϕ[2], ŷ1...) * x̂1 * rad[1]
    p5 = p2 + AngleAxis(ϕ[1], ŷ2...) * x̂2 * rad[2]
    p7 = p2 + AngleAxis(ϕ[2], ŷ2...) * x̂2 * rad[2]
    n_offset = length(eM.point)
    i1 = SVector{4,Int64}(1,3,4,6) + n_offset
    i2 = SVector{4,Int64}(3,2,5,7) + n_offset
    i3 = SVector{4,Int64}(3,4,6,7) + n_offset
    i4 = SVector{4,Int64}(4,3,5,7) + n_offset
    i5 = SVector{3,Int64}(4,6,7) + n_offset
    i6 = SVector{3,Int64}(4,7,5) + n_offset
    append!(eM.point, [p1, p2, p3, p4, p5, p6, p7])
    append!(eM.tri, [i5, i6])
    append!(eM.tet, [i1, i2, i3, i4])
    ϵ = zeros(7)
    ϵ[1:3] .= -1.0
    if is_open[1]
        ϵ[1] = 0.0
        push!(eM.tri, SVector{3,Int64}(1,6,4) + n_offset)
    end
    if is_open[2]
        ϵ[2] = 0.0
        push!(eM.tri, SVector{3,Int64}(2,5,7) + n_offset)
    end
    append!(eM.ϵ, ϵ)
    return nothing
end

function create_swept_mesh(fun_gen_2, lr, rad, num_ϕ=4, is_open::Bool=true; rot_half::Bool=true)
	l_lr = length(lr)
	if isa(rad, Float64)
		rad = zeros(l_lr) .+ rad
	else
		(l_lr == length(rad)) || error("the length of lr and length of rad must be the same")
	end
    eM = eMesh{Tri,Tet}()
    Δ_ϕ = 2 * pi / num_ϕ
    rad = rad ./ cos(Δ_ϕ / 2)
	n_θ = l_lr - 1
    for k_θ = 1:n_θ
        θ = (lr[k_θ], lr[k_θ + 1])
		rad_k = (rad[k_θ], rad[k_θ + 1])
        for k_ϕ = 1:num_ϕ
            ϕ_0 = Δ_ϕ * (k_ϕ  - 0.5 * rot_half)
            ϕ_1 = ϕ_0 + Δ_ϕ
            ϕ = (ϕ_0, ϕ_1)
            if is_open
                is_open_ = (k_θ==1, k_θ==n_θ)
            else
                is_open_ = (false, false)
            end
            add_rot_sym_segment!(eM, fun_gen_2, θ, ϕ, rad_k, is_open_)
        end
    end
	remove_degenerate!(eM)
    mesh_repair!(eM)
    return eM
end



# function output_side_points(v1::Float64, the_fun::Function, width::Float64)
#     p_d0, dir_1, dir_2 = the_fun(v1)
#
#     # check direction vectors
#     (norm(dir_1) ≈ 1.0) || error("direction 1 does not have a magnitude of 1")
#     (norm(dir_2) ≈ 1.0) || error("direction 2 does not have a magnitude of 1")
#     (norm(dot(dir_1, dir_2)) < 1.0e-14) || error("direction vectors are not orthonormal")
#
#     dir_3 = cross(dir_1, dir_2)
#     delta_y = dir_2 * width
#     delta_z = dir_3 .* width
#
#     p1 = p_d0 - delta_y + delta_z
#     p2 = p_d0 + delta_y + delta_z
#     p5 = p_d0 - delta_y - delta_z
#     p6 = p_d0 + delta_y - delta_z
#
#     return p1, p2, p5, p6
# end
#
# function create_swept_segment(width::Float64, v0::Float64, v1::Float64, the_fun::Function, open_neg_z::Bool, open_pos_z::Bool)
#     p1, p2, p3, p4 = output_side_points(v0, the_fun, width)
#     p5, p6, p7, p8 = output_side_points(v1, the_fun, width)
#     p_neg_z = (p1 + p2 + p3 + p4) ./ 4
#     p_pos_z = (p5 + p6 + p7 + p8) ./ 4
#     p_center = (p_neg_z + p_pos_z) ./ 4
#     eM_box = eMesh{Tri,Tet}()
#     append!(eM_box.ϵ, zeros(8))
#     push!(eM_box.ϵ, -1.0)
#     push!(eM_box.ϵ, -1.0)
#     push!(eM_box.ϵ, -1.0)
#     push!(eM_box.point, p1, p2, p3, p4, p5, p6, p7, p8, p_center, p_neg_z, p_pos_z)
#
#     # add non-z triangles
#     push!(eM_box.tri, SVector{3,Int64}(1,5,7))  # neg x
#     push!(eM_box.tri, SVector{3,Int64}(1,7,3))
#     push!(eM_box.tri, SVector{3,Int64}(2,4,8))  # pos x
#     push!(eM_box.tri, SVector{3,Int64}(2,8,6))
#     push!(eM_box.tri, SVector{3,Int64}(1,2,6))  # neg y
#     push!(eM_box.tri, SVector{3,Int64}(1,6,5))
#     push!(eM_box.tri, SVector{3,Int64}(3,7,8))  # pos y
#     push!(eM_box.tri, SVector{3,Int64}(3,8,4))
#
#     if open_neg_z
#         push!(eM_box.tri, SVector{3,Int64}(1,3,4))  # neg z
#         push!(eM_box.tri, SVector{3,Int64}(1,4,2))  # neg z
#     else
#         push!(eM_box.tri, SVector{3,Int64}(1,3,10))
#         push!(eM_box.tri, SVector{3,Int64}(3,4,10))
#         push!(eM_box.tri, SVector{3,Int64}(4,2,10))
#         push!(eM_box.tri, SVector{3,Int64}(2,1,10))
#     end
#
#     if open_pos_z
#         push!(eM_box.tri, SVector{3,Int64}(5,6,8))  # pos z
#         push!(eM_box.tri, SVector{3,Int64}(5,8,7))  # pos z
#     else
#         push!(eM_box.tri, SVector{3,Int64}(5,6,11))
#         push!(eM_box.tri, SVector{3,Int64}(6,8,11))
#         push!(eM_box.tri, SVector{3,Int64}(8,7,11))
#         push!(eM_box.tri, SVector{3,Int64}(7,5,11))
#     end
#
#     for k = 1:n_tri(eM_box)
#         push!(eM_box.tet, SVector{4,Int64}(9, eM_box.tri[k]...))
#     end
#
#     return eM_box
# end
#
# function make_swept_mesh_closed(width::Float64, the_fun::Function, lr)
#     eM_box = eMesh{Tri,Tet}()
#     n_segment = length(lr)
#     for k = 1:n_segment
#         eM_1 = create_swept_segment(width, lr[k], lr[mod1(k+1, n_segment)], the_fun, false, false)
#         append!(eM_box, eM_1)
#     end
#     mesh_repair!(eM_box)
#     return eM_box
# end
#
# function make_swept_mesh_open(width::Float64, the_fun::Function, lr)
#     eM_box = eMesh{Tri,Tet}()
#     n_segment = length(lr) - 1
#     for k = 1:n_segment
#         open_neg_z = (k == 1)
#         open_pos_z = (k == n_segment)
#         eM_1 = create_swept_segment(width, lr[k], lr[k+1], the_fun, open_neg_z, open_pos_z)
#         append!(eM_box, eM_1)
#     end
#     mesh_repair!(eM_box)
#     return eM_box
# end

########################

# function create_swept_segment(width::Float64, v0::Float64, v1::Float64, the_fun::Function)
#     p1, p2, p5, p6 = output_side_points(v0, the_fun, width)
#     p3, p4, p7, p8 = output_side_points(v1, the_fun, width)
#     p9 = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) ./ 8
#     eM_box = output_eMesh_box(1.0)
#     empty!(eM_box.point)
#     push!(eM_box.point, p1, p2, p3, p4, p5, p6, p7, p8, p9)
#     return eM_box
# end
#
# function make_swept_mesh(width::Float64, the_fun::Function, lr)
#     eM_box = eMesh{Tri,Tet}()
#     for k = 1:(length(lr) - 1)
#         eM_1 = create_swept_segment(width, lr[k], lr[k+1], the_fun)
#         append!(eM_box, eM_1)
#     end
#     mesh_repair!(eM_box)
#     return eM_box
# end



function f_swept_helix(θ::Float64, coil_sep::Float64)
    delta_z = 1 / (2 * pi) * coil_sep
    r = SVector{3,Float64}(cos(θ), sin(θ), θ * delta_z)
    dir_1 = normalize(SVector{3,Float64}(-sin(θ), cos(θ), delta_z))
    dir_2 = SVector{3,Float64}( cos(θ), sin(θ), 0.0)
    return r, dir_1, dir_2
end

f_swept_circle(θ::Float64) = f_swept_helix(θ, 0.0)
