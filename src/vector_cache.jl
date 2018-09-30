mutable struct VectorCache{T}
    ind_fill::Int64
    ind_max::Int64
    vec::Vector{T}
    empty::T
    VectorCache(t::T) where {T} = new{T}(-9999, 0, Vector{T}(), t)
end

function expand!(vc::VectorCache{T}) where {T}  # TODO: make this function more elegant
    ind_expand = 64
    resize!(vc.vec, vc.ind_max + ind_expand)
    for k = 1:ind_expand
        vc.vec[k + vc.ind_max] = deepcopy(vc.empty)
    end
    vc.ind_max += ind_expand
    return nothing
end

function returnNext(vc::VectorCache{T}) where {T}
    (vc.ind_fill == vc.ind_max) && expand!(vc)
    vc.ind_fill += 1
    return vc.vec[vc.ind_fill]
end

function addCacheItem!(vc::VectorCache{T}, item_T::T) where {T}
  vc.ind_fill += 1
  (vc.ind_max < vc.ind_fill) && expand!(vc)
  vc.vec[vc.ind_fill] = item_T
  return nothing
end

function Base.empty!(vc::VectorCache{T}) where {T}
  vc.ind_fill = 0
  return nothing
end

Base.isempty(vc::VectorCache{T}) where {T} = (vc.ind_fill == 0)
Base.length(vc::VectorCache{T}) where {T} = vc.ind_fill
Base.@propagate_inbounds Base.getindex(vc::VectorCache{T}, i::Int) where {T} = vc.vec[i]
