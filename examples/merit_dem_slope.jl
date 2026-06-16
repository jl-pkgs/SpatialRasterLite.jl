using SpatialRasterLite, ArchGDAL

f = "X:/rpkgs/VICTools.R/inst/database/merit/merit_China_1km_dem.tif"
prefix = "merit_China_1km"

f = "Z:/China/ALLChinaRunoff/ChinaBasins/china90_merit/merit90_china_demfill.tif"
prefix = "merit_China_90m" # 54s
# 计算坡度，输出相同 WGS84，单位 m/m
# src/terrain/dem_slope.jl

dem = SpatRaster(f)              # 读取 DEM (WGS84)
@time slope, dir = dem_slope(dem)      # slope: m/m; dir: ArcGIS D8 坡向, 几何同 dem

write_gdal(slope, "$(prefix)_slope_mm.tif")
write_gdal(dir, "$(prefix)_slope_dir.tif")
