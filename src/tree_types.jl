mutable struct bin_BB_Tree{T}
    id::Int64
    box::T
    node_1::bin_BB_Tree{T}
    node_2::bin_BB_Tree{T}
    # function bin_BB_Tree{T}() where {T <: boundingBox}
    #     return new(-9999)
    # end
    function bin_BB_Tree{T}(id::Int64, BB::T) where {T <: boundingBox}
        return new(id, BB)
    end
    # function bin_BB_Tree{T}(BB::T, node_1::bin_BB_Tree{T}, node_2::bin_BB_Tree{T}) where {T <: boundingBox}
    #     return new(-9999, BB, node_1, node_2)
    # end
    function bin_BB_Tree{AABB}(node_1::bin_BB_Tree{AABB}, node_2::bin_BB_Tree{AABB})
        # min_1, max_1 = boxMinMax(node_1.box)
        # min_2, max_2 = boxMinMax(node_2.box)
        # aabb = svSvToAABB(SVector{4,SVector{3,Float64}}(min_1, max_1, min_2, max_2))
        aabb = combineAABB(node_1.box, node_2.box)
        return new(-9999, aabb, node_1, node_2)
    end
end


mutable struct TT_Cache
    vc::vectorCache{Tuple{Int64, Int64}}
    q_a_b::Quat{Float64}
    t_a_b::SVector{3,Float64}
    function TT_Cache()
        vc = vectorCache{Tuple{Int64,Int64}}()
        return new(vc)
    end
end

isleaf(tree::bin_BB_Tree) = (tree.id == -9999)

function tree_tree_intersect(tree_1::bin_BB_Tree{T}, tree_2::bin_BB_Tree{T}, ttCache::TT_Cache) where {T <: boundingBox}
    if BB_BB_intersect(tree_1.box, tree_2.box, ttCache.q_a_b, ttCache.t_a_b)
        isleaf_1 = isleaf(tree_1)
        isleaf_2 = isleaf(tree_2)
        if isleaf_1
            if isleaf_2
                addCacheItem!(ttCache.vc, (tree_1.id, tree_2.id))
            else
                tree_tree_intersect(tree_1, tree_2.node_1, ttCache)
                tree_tree_intersect(tree_1, tree_2.node_2, ttCache)
            end
        else
            if isleaf_2
                tree_tree_intersect(tree_1.node_1, tree_2.node_1, ttCache)
                tree_tree_intersect(tree_1.node_2, tree_2.node_1, ttCache)
                tree_tree_intersect(tree_1.node_1, tree_2.node_2, ttCache)
                tree_tree_intersect(tree_1.node_2, tree_2.node_2, ttCache)
            else
                tree_tree_intersect(tree_1.node_1, tree_2, ttCache)
                tree_tree_intersect(tree_1.node_2, tree_2, ttCache)
            end
        end
    end
    return nothing
end
