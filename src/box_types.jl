abstract type boundingBox end

struct AABB <: boundingBox
    c::SVector{3,Float64}
    e::SVector{3,Float64}
    AABB(c, e) = new(c, e)
end

struct OBB <: boundingBox
    c::SVector{3,Float64}
    e::SVector{3,Float64}
    q::Quat{Float64}
    OBB(c, e, q) = new(c, e, q)
end

boxMinMax(a::AABB) = a.c - a.e, a.c + a.e
function combineAABB(a::AABB, b::AABB)
    min_1, max_1 = boxMinMax(a)
    min_2, max_2 = boxMinMax(b)
    return svSvToAABB(SVector{4,SVector{3,Float64}}(min_1, max_1, min_2, max_2))
end
