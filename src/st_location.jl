export _st_location, st_location, st_location_exact


"""
    st_location(r::Raster, points::Vector{Tuple{T,T}})

return the overlaping indexes `inds`, and corresponding (i,j)

# Examples
```julia
inds, locs = st_location(r, points)
```
"""
function _st_location(rastsize::RasterSize, x::T, y::T) where {T<:Real}
  (; b, cellsize, nx, ny) = rastsize
  (b.xmin <= x <= b.xmax && b.ymin <= y <= b.ymax) || return nothing
  # Match terra: internal ties go east/south; include the outer edges.
  i = x == b.xmax ? nx : floor(Int, (x - b.xmin) * (nx / (b.xmax - b.xmin))) + 1
  j = y == b.ymin ? ny : floor(Int, (b.ymax - y) * (ny / (b.ymax - b.ymin))) + 1
  (1 <= i <= nx && 1 <= j <= ny) || return nothing
  i, cellsize[2] > 0 ? ny - j + 1 : j
end

function st_location(rastersize::RasterSize, point::P) where {
  T<:Real,P<:Union{Tuple{T,T},Point{T}}}
  _st_location(rastersize, get_x(point), get_y(point))
end

function st_location(rastersize::RasterSize, points::Vector{P}; rm_empty::Bool=false) where {
  T<:Real,P<:Union{Tuple{T,T},Point{T}}}
  map(p -> _st_location(rastersize, get_x(p), get_y(p)), points)
  # return rm_empty ? _rm_empty(locs) : locs
end


## 全部适用于raster
function st_location(ra::AbstractSpatRaster, args...; kw...)
  rastersize = RasterSize(ra)
  st_location(rastersize, args...; kw...)
end


# - `tol`: (0.5 + tol) * cellsize
function st_location_exact(lon::AbstractVector, lat::AbstractVector, points::Vector{Tuple{T,T}};
  rm_empty::Bool=false, cellsize=nothing, tol=1e-2) where {T<:Real}
  isnothing(cellsize) && (cellsize = st_cellsize(lon, lat))
  cellx, celly = cellsize
  
  map(p -> findnear(p, lon, lat; cellx, celly, tol), points)
  # return rm_empty ? _rm_empty(locs) : locs
end

# function xy2ij(x::T, y::T, b::bbox, cellsize) where {T<:Real}
#   cellx, celly = cellsize
#   i = floor(Int, (x - b.xmin) / cellx)
#   if celly < 0
#     j = floor(Int, (b.ymax - y) / abs(celly))
#   else
#     j = floor(Int, (y - b.ymin) / celly)
#   end
#   return i, j
# end

# function xy2ij(point::Tuple{T,T}, ra::SpatRaster) where {T<:Real}
#   x, y = point
#   xy2ij(x, y, st_bbox(ra), st_cellsize(ra))
# end
