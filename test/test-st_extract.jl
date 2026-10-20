# using Ipaper, Ipaper.sf, Test, ArchGDAL

@testset "st_extract" begin
  points = [
    (110.7, 32.25),
    (110.72, 32.27),
    (110.71354166666667, 32.484375),
    (111.7, 32.25)
  ]
  f = guanshan_dem
  ra = rast(f)
  inds, vals = st_extract(ra, points)
  @test inds == 1:4
  @test length(vals) == 4
  @test isnan(vals[end])
  # r2 = st_resample(ra; fact=10)
  # @test size(r2) == (16, 12, 1)
end

@testset "resample_first" begin
  r = resample_first(rand(10, 10, 4), fact=2)
  @test size(r) == (5, 5, 4)
  @test isa(r, SubArray)

  r = resample_first(rand(10, 10, 4), fact=2, deepcopy=true)
  @test !isa(r, SubArray)  
end

@testset "st_location" begin
  lon = 0.5:1:9.5
  lat = 0.5:1:9.5
  nlon, nlat = length(lon), length(lat)
  b = st_bbox(lon, lat)
  A = rand(nlon, nlat)
  ra = SpatRaster(A, b)

  points = [
    (0., 0.),
    (0.5, 0.5),
    (0.1, 0.9),
    (1.0, 1.0), # 处于边界上
    (9.5, 9.9),
    (10, 10)    # 外边界，纳入栅格
  ]
  locs = st_location(ra, points)
  # @test which_notnull(locs) == 1:5
  @test locs == [(1, 10), (1, 10), (1, 10), (2, 10), (10, 1), (10, 1)]

  locs = st_location_exact(lon, lat, points)
  # @test which_notnull(locs) == 1:6
end
