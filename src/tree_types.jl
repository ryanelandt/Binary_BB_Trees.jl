mutable struct bin_BB_Tree{T}
    id::Int64
    box::T
    node_1::bin_BB_Tree{T}
    node_2::bin_BB_Tree{T}
    function bin_BB_Tree{T}() where {T <: boundingBox}
        return new(-9999)
    end
end
