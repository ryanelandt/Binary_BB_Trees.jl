unit_cube_points = Binary_BB_Trees.basicBoxPoints()

@testset "unit cube" begin
    @test unit_cube_points[1] == SVector{3,Float64}(-1.0, -1.0, -1.0)
    @test unit_cube_points[2] == SVector{3,Float64}(+1.0, -1.0, -1.0)
    @test unit_cube_points[3] == SVector{3,Float64}(-1.0, +1.0, -1.0)
    @test unit_cube_points[4] == SVector{3,Float64}(+1.0, +1.0, -1.0)
    @test unit_cube_points[5] == SVector{3,Float64}(-1.0, -1.0, +1.0)
    @test unit_cube_points[6] == SVector{3,Float64}(+1.0, -1.0, +1.0)
    @test unit_cube_points[7] == SVector{3,Float64}(-1.0, +1.0, +1.0)
    @test unit_cube_points[8] == SVector{3,Float64}(+1.0, +1.0, +1.0)
    five_tets = Binary_BB_Trees.cubeToFiveTets(unit_cube_points)
    @test 8.0 ≈ sum(volume.(five_tets))
end

@testset "orientation: cube face / triangle" begin
    i_box_faces = Binary_BB_Trees.outputOrientedBoxFaces()
    tup_dir = (
        SVector{3,Float64}(-1.0,  0.0,  0.0),
        SVector{3,Float64}(+1.0,  0.0,  0.0),
        SVector{3,Float64}( 0.0, -1.0,  0.0),
        SVector{3,Float64}( 0.0, +1.0,  0.0),
        SVector{3,Float64}( 0.0,  0.0, -1.0),
        SVector{3,Float64}( 0.0,  0.0, +1.0)
    )
    for (k, dir_k) = enumerate(tup_dir)
        v_tri = Vector{SVector{3,SVector{3,Float64}}}()
        Binary_BB_Trees.twoTriangles!(v_tri, unit_cube_points[i_box_faces[k]])
        @test dir_k ≈ triangleNormal(v_tri[2])
        @test dir_k ≈ triangleNormal(v_tri[1])
    end
end
