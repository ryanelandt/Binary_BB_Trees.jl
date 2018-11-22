
function top_down(point::Vector{SVector{3,Float64}}, vec_tri_tet::Vector{SVector{N,Int64}}) where {N}
    (N == 3) || (N == 4) || error("N is $N. N should be 3 for triangle mesh or 4 for tetrahedral mesh")

    sv_aabb = [svSvToAABB(point[k]) for k = vec_tri_tet]
    bb_aabb = [bin_BB_Tree(k, sv_aabb[k]) for k = 1:length(sv_aabb)]
    return recursive_top_down(bb_aabb)
end

function bin_BB_Tree_Bound(bb_aabb)
    box = bb_aabb[1].box
    for k = 2:length(bb_aabb)
        box = combineAABB(box, bb_aabb[k].box)
    end
    return box
end

function recursive_top_down(bb_aabb::Vector{bin_BB_Tree})
    n = length(bb_aabb)
    (n == 1) && (return bb_aabb[1])
    (n == 2) && (return bin_BB_Tree(bb_aabb[1], bb_aabb[2]))

    final_AABB = bin_BB_Tree_Bound(bb_aabb)
    _, mi_ = findmax(final_AABB.e)

    all_c = [bb_aabb[k].box.c[mi_] for k = 1:n]
    s_perm = sortperm(all_c)
    n_mid = Int64(ceil(n / 2))
    i_pos = s_perm[1:(n_mid - 1)]
    i_neg = s_perm[n_mid:n]

    bb_aabb_1 = bb_aabb[i_pos]
    bb_aabb_2 = bb_aabb[i_neg]
    tree_1 = recursive_top_down(bb_aabb_1)
    tree_2 = recursive_top_down(bb_aabb_2)
    return bin_BB_Tree(tree_1, tree_2)
end
