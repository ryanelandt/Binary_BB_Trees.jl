# struct SurfaceVolumeMesh
#     point::Vector{SVector{3,Float64}}
#     tri::Vector{SVector{3,Int64}}
#     tet::Vector{SVector{4,Int64}}
#     function SurfaceVolumeMesh()
#         point = Vector{SVector{3,Float64}}()
#         tri = Vector{SVector{3,Int64}}()
#         tet = Vector{SVector{4,Int64}}()
#         return new(point, tri, tet)
#     end
#     SurfaceVolumeMesh(point::Vector{SVector{3,Float64}}, tri::Vector{SVector{3,Int64}}, tet::Vector{SVector{4,Int64}}) = new(point, tri, tet)
# end
