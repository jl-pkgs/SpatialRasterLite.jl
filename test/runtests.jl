using Test, SpatialRasterLite, ArchGDAL

import Random: seed!
set_seed(seed) = seed!(seed)

# include("hydro/test-flowdir.jl")
include("test-sun_angle.jl")
include("test-intersect.jl")
include("test-distance.jl")

include("terrain/test-dem_slope.jl")

# include("Interp/test-interp.jl")
# include("Interp/test-tps.jl")
# include("Interp/test-angle.jl")
# include("Interp/test-bilinear.jl")
# include("Interp/test-bilinear_irregular.jl")

include("test-gdal_polygonize.jl")
include("test-bbox.jl")
include("test-rast.jl")
include("test-st_extract.jl")
include("test-extract-boundaries.jl")
include("test-st_mosaic.jl")
include("test-write_gdal.jl")

# println(dirname(@__FILE__))
# println(pwd())
# cd(dirname(@__FILE__)) do
