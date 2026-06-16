export dem_angle_MaxElevation
export dem_slope


# ArcGIS D8 flow-direction codes (最陡下降方向)
#   32  64  128
#   16   ·    1
#    8   4    2
"""
    dem_slope(dem::SpatRaster; R=6378.388)

由 WGS84 高程栅格计算**最陡下降方向**的坡度与坡向 (D8)。

对每个格点, 计算其指向 8 个相邻格点的下降梯度 `(z0 - z_i) / dist_i`, 取最大者:
- `slope`: 该最陡下降方向的坡度, 单位 **m/m** (坡度角的正切, rise/run); 浮点,
  边界格点为 `NaN`;
- `dir`:   最陡下降方向, 按 ArcGIS flow direction 编码 (E=1, SE=2, S=4, SW=8,
  W=16, NW=32, N=64, NE=128), 以 `UInt8` 存储; 边界、洼地/平地等无下降方向记 `0`。

经向、纬向及对角的地面格距由经纬度换算为实际距离 [m], 随纬度逐行变化。两个输出
`SpatRaster` 与 `dem` 几何范围相同。

# Arguments
- `dem`: 高程栅格, 高程单位 m, 经纬度单位为度。
- `R`: 地球半径 [km], 用于地面距离换算。

# Return
`(; slope, dir)`, 均为 `SpatRaster` (`slope` 为浮点, `dir` 为 `UInt8`)。
"""
function dem_slope(dem::SpatRaster{T,2}; R=6378.388) where {T}
  lon, lat = st_dims(dem)
  cellx, celly = abs.(st_cellsize(dem))
  FT = float(T)
  slope = fill(FT(NaN), size(dem.A))
  dir = zeros(UInt8, size(dem.A))  # 0 = 无下降方向 (边界/洼地/nodata)
  # 默认 reverse_lat=true → lat 降序 → j+1 指向南; lat 升序时南北向编码互换
  lat_desc = length(lat) < 2 || lat[2] < lat[1]
  # 距离/几何用 Float64 (与经纬度一致), 仅结果写回 FT 数组
  _dem_slope!(slope, dir, dem.A, lon, lat, Float64(cellx), Float64(celly), Float64(R), lat_desc)

  rs = SpatRaster(slope, dem); rs.name = "slope"
  rd = SpatRaster(dir, dem); rd.name = "flowdir"
  (; slope=rs, dir=rd)
end

# 函数屏障: 在此对 z 的具体类型特化, 消除 SpatRaster.A 抽象字段引入的类型不稳定
function _dem_slope!(slope::Matrix{FT}, dir::Matrix{UInt8}, z::AbstractMatrix,
  lon, lat, cellx::Float64, celly::Float64, R::Float64, lat_desc::Bool) where {FT}
  nlon, nlat = size(z)
  # 邻域 ArcGIS 编码. lat 降序 (默认): j+1=南, j-1=北; 升序时南北互换
  cS, cN = lat_desc ? (0x04, 0x40) : (0x40, 0x04)
  cSE, cNE = lat_desc ? (0x02, 0x80) : (0x80, 0x02)
  cSW, cNW = lat_desc ? (0x08, 0x20) : (0x20, 0x08)
  cE, cW = 0x01, 0x10

  p = Progress(nlat - 2)
  @inbounds @threads for j in 2:nlat-1
    next!(p)
    φ, λ = lat[j], lon[1]
    # 经向、纬向、对角的地面格距随纬度变化, 逐行计算 [km] -> [m]
    dx = earth_dist((λ, φ), (λ + cellx, φ); R) * 1000
    dy = earth_dist((λ, φ), (λ, φ + celly); R) * 1000
    dxy = sqrt(dx^2 + dy^2)

    for i in 2:nlon-1
      z0 = z[i, j]
      isnan(z0) && continue
      # 8 邻域下降梯度, 取最陡下降者及其方向 (NaN 邻域经 -Inf 初值自动跳过)
      smax = -Inf; cd = 0x00
      g = (z0 - z[i+1, j]) / dx;     g > smax && (smax = g; cd = cE)
      g = (z0 - z[i-1, j]) / dx;     g > smax && (smax = g; cd = cW)
      g = (z0 - z[i, j+1]) / dy;     g > smax && (smax = g; cd = cS)
      g = (z0 - z[i, j-1]) / dy;     g > smax && (smax = g; cd = cN)
      g = (z0 - z[i+1, j+1]) / dxy;  g > smax && (smax = g; cd = cSE)
      g = (z0 - z[i-1, j+1]) / dxy;  g > smax && (smax = g; cd = cSW)
      g = (z0 - z[i+1, j-1]) / dxy;  g > smax && (smax = g; cd = cNE)
      g = (z0 - z[i-1, j-1]) / dxy;  g > smax && (smax = g; cd = cNW)

      if smax > 0
        slope[i, j] = smax
        dir[i, j] = cd
      elseif isfinite(smax)
        slope[i, j] = FT(0)   # 洼地/平地: 有邻域但无下降, dir 保持 0
      end
      # smax = -Inf (邻域全 NaN): slope 保持 NaN, dir 保持 0
    end
  end
  return nothing
end


"slope in radian"
function dem_slope(p0::Point3{T}, p1::Point3{T}) where {T}
  dl = earth_dist(p0, p1) * 1000 # 水平面上的距离, [km] to [m]
  dz = p1.z - p0.z # [m]
  atan(dz / dl) # radians, [-90°, 90°]
end

function dem_slope(p0::Point3{T}, Points::Vector{Point3{T}}) where {T}
  map(p1 -> dem_slope(p0, p1), Points) # αs, H = pi/2 - maximum(αs)
end


## 提前算好，各个方向的最大坡度
function dem_angle_MaxElevation(elev::SpatRaster, p0::Point{T};
  δψ=3, radian=2.0, rastersize::RasterSize) where {T}

  z0 = st_extract(elev, [(p0.x, p0.y)]).value[1] # 
  P0 = Point3(p0.x, p0.y, z0)

  ψs = δψ/2:δψ:360 # 天文学方位角
  # ψs =[180.0]
  map(Φ_sun -> begin
      l = dLine(; origin=p0, azimuth=Φ_sun, length=radian) # 200km^2
      points = intersect(elev, l, rastersize)
      length(points) == 0 && return NaN

      αs = dem_slope(P0, points) # 最大仰角对应的[radian]
      maximum(αs)
    end, ψs)
end


function dem_angle_MaxElevation(elev::SpatRaster; δψ=3, radian=2.0)
  # cellsize = st_cellsize(elev)
  rastersize = RasterSize(elev)
  lon, lat = st_dims(elev)
  nlon, nlat = length(lon), length(lat)
  # nlon = 20

  ψs = δψ/2:δψ:360 # 天文学方位角
  N = length(ψs)
  R = zeros(nlon, nlat, N)
  
  p = Progress(nlon)
  @inbounds @threads for i in 1:nlon
    next!(p)
    for j in 1:nlat
      p0 = Point(lon[i], lat[j])
      R[i, j, :] .= dem_angle_MaxElevation(elev, p0; δψ, radian, rastersize)
    end
  end
  rast(R, st_bbox(elev))
end
