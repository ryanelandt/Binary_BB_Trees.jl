struct blob
    k::Int64
    n_below::Int64
    cost::Float64
    neighbor::Set{Int64}
    bin_BB_Tree::bin_BB_Tree{AABB}
    function blob(k, n_below::Int64, neighbor_in, bb_tree::bin_BB_Tree{AABB})
        cost = blobCost(bb_tree.box, n_below)
        return new(k, n_below, cost, Set{Int64}(neighbor_in), bb_tree)
    end
end

function createBlobDictionary(point::Vector{SVector{3,Float64}}, vec_tri_tet::Vector{SVector{N,Int64}}) where {N}
    vec_neighbor = extractTriTetNeighborInformation(vec_tri_tet)
    dict_blob = Dict{Int64,blob}()
    for (k, ind_k) = enumerate(vec_tri_tet)
        aabb_k = svSvToAABB(point[ind_k])
        tree_k = bin_BB_Tree{AABB}(k, aabb_k)
        dict_blob[k] = blob(k, 1, vec_neighbor[k], tree_k)
    end
    return dict_blob
end

function extractTriTetNeighborInformation(vec_tri_tet::Vector{SVector{N,Int64}}) where {N}
    dict_vert_pair = createSharedEdgeFaceDict(vec_tri_tet)
    vec_neighbor = [Vector{Int64}() for _ = vec_tri_tet]
    for key_k = keys(dict_vert_pair)
        i_pair = dict_vert_pair[key_k]
        if (N == 3) || all(i_pair .!= -9999)  # tets will not have a partner if on the boundary
            push!(vec_neighbor[i_pair[1]], i_pair[2])
            push!(vec_neighbor[i_pair[2]], i_pair[1])
        end
    end
    return vec_neighbor
end

function createSharedEdgeFaceDict(vec_tri_tet::Vector{SVector{N,Int64}}) where {N}
    dict_vert_pair = Dict{SVector{N-1, Int64}, MVector{2,Int64}}()
    for (k_tet_tri, i_vertices) = enumerate(vec_tri_tet)
        for j = 1:N
            sv_new = sortEdgeFace(i_vertices, j)
            if haskey(dict_vert_pair, sv_new)
                mv_existing = dict_vert_pair[sv_new]
                @assert(mv_existing[2] .== -9999, "three tri/tet $N share the same edge something is wrong")
                mv_existing[2] = k_tet_tri
            else
                dict_vert_pair[sv_new] = MVector{2,Int64}(k_tet_tri, -9999)
            end
        end
    end
    return dict_vert_pair
end

function blobCost(aabb::AABB, n_below::Int64)
    # TODO: make costs depend of final BB size
    cA = 1.0
    cV = 1.0
    V = 0.0
    V += n_below * log(2 * n_below)
    V += cA * boxArea(aabb)
    V += cV * boxVolume(aabb)
    return V  # PriorityQueue returns lowest first
end

function doCombineBlob(dict_blob, a::blob, b::blob, k_next::Int64)
    delete!(a.neighbor, b.k)
    delete!(b.neighbor, a.k)
    c_neighbor = union(a.neighbor, b.neighbor)  # add to neighbor_c
    tree_c = bin_BB_Tree{AABB}(a.bin_BB_Tree, b.bin_BB_Tree)
    n_below_c = a.n_below + b.n_below
    blob_c = dict_blob[k_next] = blob(k_next, n_below_c, c_neighbor, tree_c)
    return k_next + 1, blob_c
end

function calcMarginalCost(a::blob, b::blob)
    c_cost = blobCost(combineAABB(a.bin_BB_Tree.box, b.bin_BB_Tree.box), a.n_below + b.n_below)
    return c_cost - a.cost - b.cost
end

function refreshCostQueue!(dict_blob, pq_delta_cost, blob_a::blob, blob_c::blob)
    a_k = blob_a.k
    for b_k = blob_a.neighbor  # for each neighbor in blob_a
        delete!(pq_delta_cost, minmax(b_k, a_k))  # delete old cost
        pq_delta_cost[minmax(b_k, blob_c.k)] = calcMarginalCost(dict_blob[b_k], blob_c)  # add new cost
        blob_k = dict_blob[b_k]
        (a_k in blob_k.neighbor) || error("something is wrong")
        delete!(blob_k.neighbor, a_k)
        push!(blob_k.neighbor, blob_c.k)
    end
    delete!(dict_blob, a_k)
end

function createBlobPriorityQueue(dict_blob)
    pq_delta_cost = PriorityQueue{Tuple{Int64,Int64},Float64}()
    for k_a = keys(dict_blob)
        blob_a = dict_blob[k_a]
        for k_b = blob_a.neighbor
            (k_a < k_b) && (pq_delta_cost[(k_a, k_b)] = calcMarginalCost(blob_a, dict_blob[k_b]))
        end
    end
    return pq_delta_cost
end

function bottomUp!(dict_blob, pq_delta_cost)
    k_next = length(dict_blob) + 1
    while !isempty(pq_delta_cost)
        key_lowest, val_lowest = dequeue_pair!(pq_delta_cost)  # automatically removes entry
        a_k, b_k = key_lowest
        blob_a = dict_blob[a_k]
        blob_b = dict_blob[b_k]
        k_next, blob_c = doCombineBlob(dict_blob, blob_a, blob_b, k_next)
        refreshCostQueue!(dict_blob, pq_delta_cost, blob_a, blob_c)  # deal with blob_a
        refreshCostQueue!(dict_blob, pq_delta_cost, blob_b, blob_c)  # deal with blob_b
    end
    return nothing
end

function triTetMeshToTreeAABB(point::Vector{SVector{3,Float64}}, vec_tri_tet::Vector{SVector{N,Int64}}) where {N}
    dict_blob = createBlobDictionary(point, vec_tri_tet)
    pq_delta_cost = createBlobPriorityQueue(dict_blob)
    bottomUp!(dict_blob, pq_delta_cost)
    all_tree = collect(values(dict_blob))
    if (length(all_tree) != 1)
        @warn "mesh is noncontinuous"
        return all_tree  # TODO do something more with this output
    else
        return all_tree[1].bin_BB_Tree
    end
end
