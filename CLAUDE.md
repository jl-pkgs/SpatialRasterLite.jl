# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`SpatialRasterLite.jl` is a lightweight Julia package for simple spatial rasters in WGS84
projection. The core type `SpatRaster` (aliased `rast`) deliberately avoids depending on heavy
geospatial packages, so the package itself is dependency-light (`using SpatialRasterLite` loads
in ~50ms). It borrows class and function naming conventions from R's `terra` and `sf` packages
(`st_*`, `read_gdal`, `write_gdal`, etc.). It exists to back spatial interpolation in
`SpatInterp.jl` and `MGWR.jl`.

## Commands

```bash
# Run the full test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# Run a single test file (must load deps first; tests assume these are in scope)
julia --project=. -e 'using Test, SpatialRasterLite, ArchGDAL; include("test/test-bbox.jl")'

# REPL with the project active
julia --project=.

# Time the load
julia --project=. -e "@time using SpatialRasterLite"
```

CI (`.github/workflows/CI.yml`) runs on Ubuntu against Julia 1.11 and 1.12 via
`julia-actions/julia-runtest`.

## Architecture

### Core data type and the GDAL boundary

The package is split across two halves that communicate through **function stubs declared but not
implemented in the main module**:

- `src/SpatialRasterLite.jl` declares empty functions (`read_gdal`, `write_gdal`, `gdalinfo`,
  `gdal_nodata`, `bandnames`, `gdal_polygonize`, `nband`, …) and exports them. The main module
  does NOT depend on ArchGDAL.
- `ext/SpatArchGDALExt/` is a **package extension** (weakdep on `ArchGDAL`) that provides the real
  implementations by `import SpatialRasterLite: read_gdal, write_gdal, …` and adding methods.

This is the key structural fact: **any actual file I/O (`read_gdal`, `write_gdal`, `gdalinfo`,
polygonize) only works once `ArchGDAL` is loaded.** Tests and downstream code do `using ArchGDAL`
to trigger the extension. When editing I/O behavior, the stub/signature lives in `src/IO.jl` but
the implementation lives in `ext/SpatArchGDALExt/IO.jl`. A second extension hook exists for
`GeoDataFrames` (`read_sf`/`write_sf` stubs).

### `SpatRaster{T,N}` (src/SpatRaster.jl)

The central mutable struct holds the data array `A`, a `bbox`, `cellsize`, explicit `lon`/`lat`
coordinate vectors, plus optional `time`, `bands`, `name`, `nodata`. Conventions to remember:

- **Latitude is stored in descending order by default** (`reverse_lat=true`); `cellsize[2]` is
  negative to signal this (see `bbox2cellsize` returning `-celly`).
- Outer constructors derive `lon`/`lat`/`cellsize` from a `bbox` + array size via `bbox2dims` /
  `bbox2cellsize` (src/bbox.jl). There is also `SpatRaster(A, ::SpatRaster)` to rebuild with the
  same geometry, and `SpatRaster(f::String)` to read a file directly.
- Arithmetic and comparison operators (`+ - * / > < & |` …) are metaprogrammed over
  `SpatRaster` in a loop, broadcasting elementwise and returning a new raster with the same
  geometry. `getindex` returns a cropped sub-`SpatRaster` (carrying geometry) for ranges, or a
  scalar for integer indices.

### `bbox` (src/bbox.jl)

Immutable `bbox(xmin, ymin, xmax, ymax)` with a family of `bbox2*` converters. `bbox2dims` and
`bbox_overlap` are the workhorses translating between geographic extent, cellsize, and
array/coordinate indices — most cropping/mosaic/extract logic flows through them.

### Module include order (src/SpatialRasterLite.jl)

Foundational (`datatype.jl`, `bbox.jl`, `tools_Ipaper.jl`) → core (`SpatRaster.jl`, `Ops.jl`) →
`st_*` spatial operations (`st_bbox`, `st_dims`, `st_extract`, `st_location`, `st_resample`,
`st_mosaic`, `st_crop`) → `IO.jl` → `methods/intersect.jl` → `terrain/terrain.jl` (which pulls in
slope, sun angle, sun shade, sky view factor). The `hydro/` directory (D8/flow direction) exists
but is **commented out** of the module and its test — treat it as inactive.

### Naming idioms

- `st_*` = spatial operations on rasters (terra/sf style). `rast` is the public alias for
  `SpatRaster`; `nlyr = nband`, `st_read = read_gdal`, `st_write = write_gdal`.
- Bundled test data is exposed as consts: `guanshan_dem`, `guanshan_flowdir_cpp`,
  `guanshan_flowdir_gis` (GeoTIFFs under `data/`).
