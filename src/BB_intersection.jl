if !@isdefined fB
    const fA = CartesianFrame3D("A")
    const fa = CartesianFrame3D("a")
    const fb = CartesianFrame3D("b")
    const fB = CartesianFrame3D("B")
end

function BB_BB_intersect(a::AABB, b::AABB, q_a_b::Quat{Float64}, t_a_b::SVector{3,Float64})
    return BB_BB_intersect(a, b, RotMatrix{3, Float64}(q_a_b), t_a_b)
end

function BB_BB_intersect(a::AABB, b::AABB, R_a_b::RotMatrix{3, Float64}, t_a_b::SVector{3,Float64})
    t = R_a_b * b.c + t_a_b - a.c
    return BB_BB_intersect(a.e, b.e, R_a_b, t)
end

function BB_BB_intersect(a::OBB, b::OBB, q_a_b::Quat{Float64}, t_a_b::SVector{3,Float64})
    x_a_A = Transform3D(fA, fa, a.q, a.c)
    x_a_b = Transform3D(fb, fa, q_a_b, t_a_b)
    x_b_B = Transform3D(fB, fb, b.q, b.c)
    x_A_b = inv(x_a_A) * x_a_b
    x_A_B = x_A_b * x_b_B
    R = rotation(x_A_B)
    t = translation(x_A_B)
    return BB_BB_intersect(a.e, b.e, R, t)
end

s221(r::SVector{3, Float64}) = SVector{3, Float64}(r[3], r[3], r[2])
s100(r::SVector{3, Float64}) = SVector{3, Float64}(r[2], r[1], r[1])

function BB_BB_intersect(e_a::SVector{3, Float64}, e_b::SVector{3, Float64}, R_mat::RotMatrix{3, Float64}, t::SVector{3, Float64})
  # This implementation of the OBB_OBB_SAT test is based off of the test described
  # in table 4.1 on page 103 of Ericson's Real-Time Collision Detection.

  R = SMatrix{3,3,Float64,9}(R_mat)
  abs_R = abs.(R) + 1.0e-14

  R0_012 = SVector{3, Float64}(R[1], R[4], R[7])
  R1_012 = SVector{3, Float64}(R[2], R[5], R[8])
  R2_012 = SVector{3, Float64}(R[3], R[6], R[9])
  abs_R0_012 = SVector{3, Float64}(abs_R[1], abs_R[4], abs_R[7])
  abs_R1_012 = SVector{3, Float64}(abs_R[2], abs_R[5], abs_R[8])
  abs_R2_012 = SVector{3, Float64}(abs_R[3], abs_R[6], abs_R[9])

  # face test 1/2
  T_dot_L = abs.(t)
  r_a = e_a
  r_b = abs_R * e_b
  any((r_a + r_b) .< T_dot_L) && (return false)

  # face test 2/2
  T_dot_L = abs.(transpose(R) * t)  # these are transposes
  r_a = transpose(abs_R) * e_a
  r_b = e_b
  any((r_a + r_b) .< T_dot_L) && (return false)

  eb_100     = s100(e_b)
  eb_221     = s221(e_b)
  abs_R0_221 = s221(abs_R0_012)
  abs_R0_100 = s100(abs_R0_012)
  abs_R1_221 = s221(abs_R1_012)
  abs_R1_100 = s100(abs_R1_012)
  abs_R2_221 = s221(abs_R2_012)
  abs_R2_100 = s100(abs_R2_012)

  t0 = t[1]
  t1 = t[2]
  t2 = t[3]
  ea0 = e_a[1]
  ea1 = e_a[2]
  ea2 = e_a[3]

  # cross test 1/3
  T_dot_L = abs.(t2      * R1_012     - t1      * R2_012)
  r_a     =      ea1     * abs_R2_012 + ea2     * abs_R1_012
  r_b     =      eb_100 .* abs_R0_221 + eb_221 .* abs_R0_100
  any((r_a + r_b) .< T_dot_L) && (return false)

  # cross test 2/3
  T_dot_L = abs.(t0  * R2_012     -     t2  * R0_012)
  r_a     =      ea0 * abs_R2_012 +    ea2  * abs_R0_012
  r_b     =  eb_100 .* abs_R1_221 + eb_221 .* abs_R1_100
  any((r_a + r_b) .< T_dot_L) && (return false)

  # cross test 3/3
  T_dot_L =  abs.(t1 * R0_012     -     t0  * R1_012)
  r_a     =      ea0 * abs_R1_012 +    ea1  * abs_R0_012
  r_b     =  eb_100 .* abs_R2_221 + eb_221 .* abs_R2_100
  any((r_a + r_b) .< T_dot_L) && (return false)
  return true
end