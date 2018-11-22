
struct AABB
    c::SVector{3,Float64}
    e::SVector{3,Float64}
end

boxArea(a::AABB) = 2 * dot(a.e, SVector{3,Float64}(a.e[2], a.e[3], a.e[1]))
boxVolume(a::AABB) = prod(a.e)

function combineAABB(a::AABB, b::AABB)
    # TODO: put in utility

    min_1, max_1 = calc_min_max(a)
    min_2, max_2 = calc_min_max(b)
    return svSvToAABB(SVector{4,SVector{3,Float64}}(min_1, max_1, min_2, max_2))
end
