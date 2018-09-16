mutable struct bin_BB_Tree{T}
    id::Int64
    box::T
    node_1::bin_BB_Tree{T}
    node_2::bin_BB_Tree{T}
    function bin_BB_Tree{T}(id::Int64, BB::T) where {T <: boundingBox}
        return new(id, BB)
    end
    function bin_BB_Tree{AABB}(node_1::bin_BB_Tree{AABB}, node_2::bin_BB_Tree{AABB})
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

is_leaf(tree::bin_BB_Tree)     = (tree.id != -9999)
is_not_leaf(tree::bin_BB_Tree) = (tree.id == -9999)

treeDepth(t::bin_BB_Tree{AABB}) = leafNumberDepth(t, 0, 0)[1]
leafNumber(t::bin_BB_Tree{AABB}) = leafNumberDepth(t, 0, 0)[2]
function leafNumberDepth(t::bin_BB_Tree{AABB}, k_depth::Int64, k_leaf::Int64)
    if is_leaf(t)
        k_leaf = 1
    else is_not_leaf(t)
        k_depth_1, k_leaf_1 = leafNumberDepth(t.node_1, k_depth + 1, k_leaf)
        k_depth_2, k_leaf_2 = leafNumberDepth(t.node_2, k_depth + 1, k_leaf)
        k_depth = max(k_depth_1, k_depth_2)
        k_leaf = k_leaf_1 + k_leaf_2
    end
    return k_depth, k_leaf
end

function tree_tree_intersect(tree_1::bin_BB_Tree{T}, tree_2::bin_BB_Tree{T}, ttCache::TT_Cache) where {T <: boundingBox}
    if BB_BB_intersect(tree_1.box, tree_2.box, ttCache.q_a_b, ttCache.t_a_b)
        is_leaf_1 = is_leaf(tree_1)
        is_leaf_2 = is_leaf(tree_2)
        if is_leaf_1
            if is_leaf_2
                addCacheItem!(ttCache.vc, (tree_1.id, tree_2.id))
            else
                tree_tree_intersect(tree_1, tree_2.node_1, ttCache)
                tree_tree_intersect(tree_1, tree_2.node_2, ttCache)
            end
        else
            if is_leaf_2
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
