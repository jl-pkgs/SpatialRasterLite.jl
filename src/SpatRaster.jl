abstract type AbstractSpatRaster{T,N} end

"""
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SpatRaster{T,N} <: AbstractSpatRaster{T,N}
  A::AbstractArray{T,N}
  b::bbox = bbox(-180.0, -90.0, 180.0, 90.0)
  cellsize::NTuple{2,Real}
  lon::AbstractVector{<:Real}
  lat::AbstractVector{<:Real}
  time::Union{AbstractVector,Nothing} = nothing
  bands::Union{AbstractVector{String},Nothing} = nothing
  name::String = "Raster"
  nodata::Union{AbstractVector{T},Nothing} = nothing
end

force_vec(x::Union{Colon,Nothing,AbstractVector}) = x
force_vec(x) = [x]

function nband(ra::SpatRaster{T,N}) where {T,N}
  n = 1
  N == 3 && (n == size(ra.A)[end])
  N > 3 && (n = size(ra.A)[end])
  return n
end

"""
    SpatRaster(A, b::bbox; reverse_lat=true, kw...)

- `kw`: other parameters: `time`, `name`
"""
function SpatRaster(A::AbstractArray{T,N}, b::bbox=bbox(-180.0, -90.0, 180.0, 90.0);
  reverse_lat=true, time=nothing, bands=nothing, name="Raster", nodata=nothing, kw...) where {T,N}

  nlyr = N >= 3 ? size(A)[end] : 1
  if N == 3 && size(A, 3) == 1
    A = A[:, :, 1]
  end

  isscalar(nodata) && (nodata = fill(nodata, nlyr))

  cellsize = bbox2cellsize(b, size(A))
  lon, lat = bbox2dims(b; cellsize, reverse_lat)
  SpatRaster(; A, b, cellsize, lon, lat, time, bands, name, nodata, kw...)
end

function SpatRaster(A::AbstractArray, r::SpatRaster; reverse_lat=true)
  (; b, time, bands, name) = r
  if size(A)[1:2] != size(r.A)[1:2]
    lon, lat = bbox2dims(b; size=size(A), reverse_lat)
    cellsize = bbox2cellsize(b, size(A))
  else
    (; lon, lat, cellsize) = r
  end
  SpatRaster(; A, b, cellsize, lon, lat, time, bands, name) # rebuild
end

function SpatRaster(f::String; FT::Union{DataType,Nothing}=nothing, kw...)
  A = read_gdal(f)
  bands = bandnames(f)
  nodata = gdal_nodata(f)

  if !isnothing(FT) 
    A = FT.(A)
    nodata = FT.(nodata)
  end
  SpatRaster(A, st_bbox(f); bands, nodata, kw...)
end


function make_rast(; b::bbox=bbox(70, 15, 140, 55), cellsize=0.5, FT=Float64)
  lon, lat = bbox2dims(b; cellsize)
  nlon, nlat = length(lon), length(lat)
  rast(zeros(FT, nlon, nlat), b)
end


function Base.getindex(ra::AbstractSpatRaster, i, j, args...; deep=true)
  (; A, cellsize, lon, lat, time, bands, name, nodata) = ra

  cols = repeat([:], max(ndims(ra) - 2 - length(args), 0))
  inds = (i, j, args..., cols...)

  if length(args) > 0
    k = force_vec(args[1])
    !isnothing(bands) && (bands = bands[k])
    !isnothing(time) && (time = time[k])
  end

  lon, lat = st_dims(ra)
  _lon, _lat = lon[i], lat[j]
  _b = st_bbox(_lon, _lat)
  _A = deep ? A[inds...] : @view A[inds...]

  rast(; A=_A, b=_b, lon=_lon, lat=_lat, time, bands,
    cellsize, name, nodata)
end

function Base.getindex(ra::AbstractSpatRaster, i::Int, j::Int, args...)
  cols = repeat([:], max(ndims(ra) - 2 - length(args), 0))
  inds = (i, j, args..., cols...)
  ra.A[inds...]
end


Base.ndims(ra::AbstractSpatRaster) = ndims(ra.A)
Base.size(ra::AbstractSpatRaster) = size(ra.A)
Base.size(ra::AbstractSpatRaster{T,2}) where {T} = (size(ra.A)..., 1)
# Base.parent(ra::AbstractSpatRaster) = ra.A
# Base.iterate(ra::AbstractSpatRaster) = iterate(ra.A)
# Base.length(ra::AbstractSpatRaster) = length(ra.A)
# Base.size(ra::AbstractSpatRaster) = size(ra.A)
# Base.eltype(::Type{AbstractSpatRaster{T}}) where {T} = T
# Base.map(f, ra::AbstractSpatRaster) = SpatRaster(map(f, ra.A), ra)

# !note about NaN values
Base_ops = ((:Base, :+), (:Base, :-), (:Base, :*), (:Base, :/),
  (:Base, :>), (:Base, :<), (:Base, :>=), (:Base, :<=),
  (:Base, :!=),
  (:Base, :&), (:Base, :|))

for (m, f) in Base_ops
  # _f = Symbol(m, ".:", f)
  @eval begin
    $m.$f(a::AbstractSpatRaster, b::AbstractSpatRaster) = begin
      size(a) != size(b) && throw(DimensionMismatch("size mismatch"))
      SpatRaster($m.$f.(a.A, b.A), a)
    end

    $m.$f(a::AbstractSpatRaster, b::Real) = SpatRaster($m.$f.(a.A, b), a)
    $m.$f(a::Real, b::AbstractSpatRaster) = SpatRaster($m.$f.(a, b.A), b)
  end
end


import Base: ==
function ==(x::SpatRaster, y::SpatRaster)
  x.b == y.b && x.A == y.A && x.nodata == y.nodata
end


function Base.show(io::IO, x::SpatRaster)
  T = eltype(x.A)
  printstyled(io, "SpatRaster{$T}: ", color=:blue)
  printstyled(io, "$(x.name)\n", color=:green, underline=true)

  print(io, "  A        : ")
  obj_size(io, x.A)

  println(io, "  b        : $(x.b)")
  println(io, "  cellsize : $(x.cellsize)")
  println(io, "  lon, lat : $(x.lon), $(x.lat)")
  if !isnothing(x.time)
    time_beg = x.time[1]
    time_end = x.time[end]
    println(io, "  time     : $time_beg ~ $time_end, ntime=$(length(x.time))")
  else
    println(io, "  time     : $(x.time)")
  end
  println(io, "  bands    : $(x.bands)")
  print(io, "  nodata   : $(x.nodata)")
  return nothing
end



@with_kw mutable struct RasterSize
  b::bbox = bbox(-180.0, -90.0, 180.0, 90.0)
  cellsize::NTuple{2,Float64} = (1.0, 1.0)
  nx::Int = guess_nx(b, cellsize)
  ny::Int = guess_nx(b, cellsize)
end
# lon::AbstractVector{<:Real}
# lat::AbstractVector{<:Real}
function guess_nx(b::bbox, cellsize)
  cellx = cellsize[1]
  length(b.ymin+cellx/2:cellx:b.ymax)
end

function guess_ny(b::bbox, cellsize)
  celly = abs(cellsize[2])
  length(b.ymin+celly+2:celly:b.ymax)
end


function RasterSize(ra::SpatRaster)
  b = st_bbox(ra)
  cellsize = ra.cellsize
  lon, lat = st_dims(ra)
  nx, ny = length(lon), length(lat)
  RasterSize(; b, cellsize, nx, ny)
end

function RasterSize(lon::AbstractVector{T}, lat::AbstractVector{T}) where {T<:Real}
  b = st_bbox(lon, lat)
  nx, ny = length(lon), length(lat)
  cellsize = st_cellsize(lon, lat)
  RasterSize(; b, cellsize, nx, ny)
end

function st_dims(rastersize::RasterSize)
  (; b, cellsize) = rastersize
  return bbox2dims(b; cellsize)
end


rast = SpatRaster

export AbstractSpatRaster, SpatRaster, rast
export make_rast
export RasterSize
