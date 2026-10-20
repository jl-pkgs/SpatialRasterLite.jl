export st_extract

function getvalue(A::AbstractArray, i::Int, j::Int)
  cols = repeat([:], ndims(A) - 2)
  A[i, j, cols...]
end

getvalue(A::AbstractArray, ij::Tuple) = getvalue(A, ij...)
getvalue(A::AbstractArray, ::Nothing) = ndims(A) == 2 ? NaN : fill(NaN, size(A)[3:end])

## TODO: `_rm_empty`让情况变得复杂
function _rm_empty(x::Vector)
  inds = findall(!isnothing, x)
  inds, x[inds]
end

# function _is_empty(locs::Vector)
#   isempty = map(p -> (p[1] != -1 && p[2] != -1), locs) # 有一个为-1就是空
#   isempty, locs
# end

function st_extract(ra::AbstractSpatRaster{FT}, point::P;) where {
  FT,T,P<:Union{Tuple{T,T}, Point{T}}}

  getvalue(ra.A, st_location(RasterSize(ra), point))
end


# Preserve input order, including out-of-bounds points (NaN).
function st_extract(ra::AbstractSpatRaster, points::Vector{P}; combine=hcat) where {
  T<:Real,P<:Union{Tuple{T,T}, Point{T}}}

  locs = st_location(ra, points)
  lst = [getvalue(ra.A, ij) for ij in locs]
  (; index=collect(eachindex(points)), value=combine(lst...))
end
