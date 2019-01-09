
abstract type BoundingBox end

struct AABB <: BoundingBox
    c::SVector{3,Float64}
    e::SVector{3,Float64}
    function AABB(c::SVector{3,Float64}, e::SVector{3,Float64})
        return new(c, e)
    end
    function AABB(aabb::AABB)
        return new(aabb.c, aabb.e)
    end
end

struct OBB <: BoundingBox
    c::SVector{3,Float64}
    e::SVector{3,Float64}
    R::SMatrix{3,3,Float64,9}
    function OBB(c::SVector{3,Float64}, e::SVector{3,Float64}, R::SMatrix{3,3,Float64,9})
        return new(c, e, R)
    end
    function OBB(aabb::AABB)
        return new(aabb.c, aabb.e, one(SMatrix{3,3,Float64,9}))
    end
end

boxArea(a::AABB) = 2 * dot(a.e, SVector{3,Float64}(a.e[2], a.e[3], a.e[1]))
boxVolume(a::AABB) = prod(a.e)
