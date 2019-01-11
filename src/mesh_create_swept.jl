
function output_side_points(v1::Float64, the_fun::Function, width::Float64)
    p_d0, dir_1, dir_2 = the_fun(v1)

    # check direction vectors
    (norm(dir_1) ≈ 1.0) || error("direction 1 does not have a magnitude of 1")
    (norm(dir_2) ≈ 1.0) || error("direction 2 does not have a magnitude of 1")
    (norm(dot(dir_1, dir_2)) < 1.0e-14) || error("direction vectors are not orthonormal")

    dir_3 = cross(dir_1, dir_2)
    delta_y = dir_2 * width
    delta_z = dir_3 .* width

    p1 = p_d0 - delta_y + delta_z
    p2 = p_d0 + delta_y + delta_z
    p5 = p_d0 - delta_y - delta_z
    p6 = p_d0 + delta_y - delta_z

    return p1, p2, p5, p6
end

function create_swept_segment(width::Float64, v0::Float64, v1::Float64, the_fun::Function)
    p1, p2, p5, p6 = output_side_points(v0, the_fun, width)
    p3, p4, p7, p8 = output_side_points(v1, the_fun, width)
    p9 = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) ./ 8
    eM_box = output_eMesh_box(1.0)
    empty!(eM_box.point)
    push!(eM_box.point, p1, p2, p3, p4, p5, p6, p7, p8, p9)
    return eM_box
end

function make_swept_mesh(width::Float64, the_fun::Function, lr)
    eM_box = eMesh{Tri,Tet}()
    for k = 1:(length(lr) - 1)
        eM_1 = create_swept_segment(width, lr[k], lr[k+1], the_fun)
        append!(eM_box, eM_1)
    end
    mesh_repair!(eM_box)
    return eM_box
end



function f_swept_helix(θ::Float64, coil_sep::Float64)
    delta_z = 1 / (2 * pi) * coil_sep
    r = SVector{3,Float64}(cos(θ), sin(θ), θ * delta_z)
    dir_1 = normalize(SVector{3,Float64}(-sin(θ), cos(θ), delta_z))
    dir_2 = SVector{3,Float64}( cos(θ), sin(θ), 0.0)
    return r, dir_1, dir_2
end

f_swept_circle(θ::Float64) = f_swept_helix(θ, 0.0)
