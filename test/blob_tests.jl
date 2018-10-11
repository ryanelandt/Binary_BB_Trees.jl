
@testset "triTetMeshToTreeAABB" begin
    # non-closed mesh
    @test_throws ErrorException triTetMeshToTreeAABB([rand(SVector{3,Float64}) for _ = 1:3], [SVector{3,Int64}(1,2,3)])
end
