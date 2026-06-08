"""
    read_gdal(file::AbstractString, indices=1, options...)

Read a raster with ArchGDAL, decoding each band by its metadata:
`out = scale * raw + offset`, with pixels equal to the band's nodata value set to
`NaN`. `indices` selects the band(s) (an `Integer` gives a 2D array, a
range/vector a 3D array); `options...` are forwarded to `ArchGDAL.read` as a
`rows, cols` window. Both match the dimensionality of the raw reader.

The element type is **not** forced to `Float32`: a band without scale/offset
metadata keeps its raw type (e.g. a `UInt8` mask stays `UInt8`, an integer
nodata keeps its sentinel since it cannot hold `NaN`), and a scaled band keeps
its own float precision (`Float32` stays `Float32`, integers promote to their
natural float type). Works for any GDAL format exposing scale/offset/nodata,
e.g. HDF4, HDF5, netCDF and GeoTIFF.
"""
function read_gdal(file::AbstractString, indices=1:nband(file), options...)
  ArchGDAL.read(file) do ds
    A = ArchGDAL.read(ds, indices, options...)
    if indices isa Integer
      decode_band(A, ArchGDAL.getband(ds, indices))           # 2D, single band
    else
      slices = map(enumerate(indices)) do (k, b)
        decode_band(@view(A[:, :, k]), ArchGDAL.getband(ds, b))
      end
      cat(slices...; dims=3)                                   # 3D, one slice per band
    end
  end
end

# Decode one 2D band: `out = scale * raw + offset`, with nodata pixels set to NaN.
# The raw element type is kept when no scaling is needed; scaling promotes
# integers to their natural float type and keeps existing float precision.
function decode_band(raw::AbstractArray, band)
  scale = ArchGDAL.getscale(band)
  offset = ArchGDAL.getoffset(band)
  nodata = ArchGDAL.getnodatavalue(band)

  out = raw
  # scale and offset
  if (scale !== nothing && scale != 1) || (offset !== nothing && offset != 0)
    T = eltype(raw) <: AbstractFloat ? eltype(raw) : float(eltype(raw))
    s, o = T(something(scale, 1)), T(something(offset, 0))
    out = @. T(raw) * s + o
  end

  # nodata to NaN
  if nodata !== nothing && eltype(out) <: AbstractFloat
    mask = raw .== convert(eltype(raw), nodata)
    out[mask] .= convert(eltype(out), NaN)
  end
  out
end


## This part is borrowed from the GeoArrays.jl package.
# MIT License, Copyright (c) 2018 Maarten Pronk
# <https://github.com/evetion/GeoArrays.jl/blob/master/src/io.jl>

# using GeoFormatTypes, ArchGDAL
# WGS84 = convert(WellKnownText, EPSG(4326))
# GFT.val(ga.crs)
const WGS84 = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]"

# const OPTIONS_DEFAULT_TIFF = Dict(
#   "TILED" => "YES", # not work
#   "COMPRESS" => "DEFLATE"
# )
import ArchGDAL: OF_UPDATE

function gdal_setproj!(f::AbstractString, transform::Vector{Cdouble})
  # ArchGDAL.open(f, "r+") do ds
  ArchGDAL.read(f; flags=OF_UPDATE) do ds
    ## Set geotransform and crs
    ArchGDAL.GDAL.gdalsetgeotransform(ds.ptr, transform)
    ArchGDAL.GDAL.gdalsetprojection(ds.ptr, WGS84)
  end
  return nothing
end


function write_gdal(data::AbstractArray, f::AbstractString;
  nodata=nothing, options=String[], NUM_THREADS=4, BIGTIFF=false)

  dtype = eltype(data)
  shortname = find_shortname(f)
  driver = ArchGDAL.getdriver(shortname)

  width, height = size(data)[1:2]
  ndims(data) == 2 && (nbands = 1)
  ndims(data) == 3 && (nbands = size(data, 3))

  if !isnothing(nodata) && !isa(nodata, Vector)
    nodata = fill(nodata, nbands)
  end

  if (shortname == "GTiff")
    options = [options..., "COMPRESS=DEFLATE", "TILED=YES", "NUM_THREADS=$NUM_THREADS"]
    BIGTIFF && (push!(options, "BIGTIFF=YES"))
  end

  try
    convert(ArchGDAL.GDALDataType, dtype)
  catch
    dtype, data = cast_to_gdal(data)
  end

  ArchGDAL.create(f; driver, width, height, nbands, dtype, options) do dataset
    for i = 1:nbands
      band = ArchGDAL.getband(dataset, i)
      ArchGDAL.write!(band, data[:, :, i])
      !isnothing(nodata) && ArchGDAL.GDAL.gdalsetrasternodatavalue(band.ptr, nodata[i])
    end
  end
end

# only support WGS84 proj
function write_gdal(ra::AbstractSpatRaster, f::AbstractString;
  nodata=nothing, options=String[], NUM_THREADS=4, BIGTIFF=true)

  isnothing(nodata) && (nodata = ra.nodata)
  write_gdal(ra.A, f; nodata, options, NUM_THREADS, BIGTIFF)
  gdal_setproj!(f, getgeotransform(ra))

  !isnothing(ra.bands) && set_bandnames(f, ra.bands)
  return f
end

# # Slice data and replace missing by nodata
# if isa(dtype, Union) && dtype.a == Missing
#   dtype = dtype.b
#   try
#     convert(ArchGDAL.GDALDataType, dtype)
#   catch
#     dtype, data = cast_to_gdal(data)
#   end
#   nodata === nothing && (nodata = typemax(dtype))
#   m = ismissing.(data)
#   data[m] .= nodata
#   data = Array{dtype}(data)
#   use_nodata = true
# end
