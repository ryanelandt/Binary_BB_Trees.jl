
abstract type BoundingBox end

struct AABB <: BoundingBox
    c::SVector{3,Float64}
    e::SVector{3,Float64}
    AABB(c::SVector{3,Float64}, e::SVector{3,Float64}) = new(c, e)
    AABB(aabb::AABB) = new(aabb.c, aabb.e)
end

struct OBB <: BoundingBox
    c::SVector{3,Float64}
    e::SVector{3,Float64}
    R::SMatrix{3,3,Float64,9}
    OBB(c::SVector{3,Float64}, e::SVector{3,Float64}, R::SMatrix{3,3,Float64,9}) = new(c, e, R)
    OBB(aabb::AABB) = new(aabb.c, aabb.e, one(SMatrix{3,3,Float64,9}))
end

boxArea(a::AABB) = 2 * dot(a.e, SVector{3,Float64}(a.e[2], a.e[3], a.e[1]))
boxVolume(a::AABB) = prod(a.e)
