
abstract type BoundingBox end

# struct AABB <: BoundingBox
#     c::SVector{3,Float64}
#     e::SVector{3,Float64}
#     AABB(c::SVector{3,Float64}, e::SVector{3,Float64}) = new(c, e)
#     AABB(aabb::AABB) = new(aabb.c, aabb.e)
# end

struct OBB <: BoundingBox
    c::SVector{3,Float64}
    e::SVector{3,Float64}
    R::SMatrix{3,3,Float64,9}
    OBB(c::SVector{3,Float64}, e::SVector{3,Float64}, R::SMatrix{3,3,Float64,9}) = new(c, e, R)
    # OBB(c::SVector{3,Float64}, e::SVector{3,Float64}) = new(c, e, one(SMatrix{3,3,Float64,9}))
    # OBB(aabb::AABB) = new(aabb.c, aabb.e, one(SMatrix{3,3,Float64,9}))
end

# boxArea(a::BB)   where {BB <: BoundingBox} = 8 * dot(a.e, SVector{3,Float64}(a.e[2], a.e[3], a.e[1]))
# boxVolume(a::BB) where {BB <: BoundingBox} = 8 * prod(a.e)

function (::Type{BB_Type})(a::BB_Type, b::BB_Type) where {BB_Type <: BoundingBox}
    min_1, max_1 = calc_min_max(a)
    min_2, max_2 = calc_min_max(b)
    return calc_obb(SVector{4,SVector{3,Float64}}(min_1, max_1, min_2, max_2))
end
