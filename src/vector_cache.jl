mutable struct vectorCache{T}
    ind_fill::Int64
    ind_max::Int64
    vec::Vector{T}
    empty::T
    function vectorCache(t::T) where {T}
        # v = Vector{T}(undef, 1)
        # if !isassigned(v, 1)
        #     try
        #         T()
        #     catch
        #         error("User supplied type ($T) does not have a trivial constructor.")
        #     end
        # end
        return new{T}(-9999, 0, Vector{T}(), t)
    end
end

function expand!(vc::vectorCache{T}) where {T}  # TODO: make this function more elegant
    ind_expand = 64
    resize!(vc.vec, vc.ind_max + ind_expand)
    # if !isassigned(vc.vec, vc.ind_max + 1)
    for k = 1:ind_expand
        vc.vec[k + vc.ind_max] = deepcopy(vc.empty)
    end
    # end
    vc.ind_max += ind_expand
    return nothing
end

function returnNext(vc::vectorCache{T}) where {T}
    (vc.ind_fill == vc.ind_max) && expand!(vc)
    vc.ind_fill += 1
    return vc.vec[vc.ind_fill]
end

function addCacheItem!(vc::vectorCache{T}, item_T::T) where {T}
  vc.ind_fill += 1
  (vc.ind_max < vc.ind_fill) && expand!(vc)
  vc.vec[vc.ind_fill] = item_T
  return nothing
end

function Base.empty!(vc::vectorCache{T}) where {T}
  vc.ind_fill = 0
  return nothing
end

Base.isempty(vc::vectorCache{T}) where {T} = (vc.ind_fill == 0)
Base.length(vc::vectorCache{T}) where {T} = vc.ind_fill
Base.@propagate_inbounds Base.getindex(vc::vectorCache{T}, i::Int) where {T} = vc.vec[i]
