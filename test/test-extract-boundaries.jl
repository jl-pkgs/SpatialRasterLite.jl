@testset "terra extraction boundaries" begin
  ra = rast(reshape(1:9, 3, 3), bbox(0, 0, 3, 3))
  points = [(0.5, 2.5), (1.0, 1.5), (1.5, 1.0), (1.0, 1.0),
            (1.5, 3.0), (3.0, 1.5), (1.5, 0.0), (0.0, 1.5),
            (0.0, 3.0), (3.0, 3.0), (0.0, 0.0), (3.0, 0.0),
            (-0.1, 1.0), (3.1, 1.0), (1.0, -0.1), (1.0, 3.1),
            (NaN, 1.0), (Inf, 1.0), (1.0, -Inf)]
  # terra::extract(..., method="simple"), ordered by input point.
  expected = [1, 5, 8, 8, 2, 6, 8, 4, 1, 3, 7, 9, fill(NaN, 7)...]
  result = st_extract(ra, points)
  @test result.index == eachindex(points)
  @test isequal(vec(result.value), expected)
  @test isequal([st_extract(ra, p) for p in points], expected)
  @test st_extract(ra, Point(1.0, 1.0)) == 8
  @test isequal(st_extract(ra, Point.(first.(points), last.(points))).value, result.value)
  @test isempty(st_extract(ra, Tuple{Float64,Float64}[]).value)

  layers = rast(cat(ra.A, ra.A .+ 100; dims=3), ra.b)
  expected_layers = permutedims(hcat(expected, expected .+ 100))
  @test isequal(st_extract(layers, points).value, expected_layers)
  @test st_extract(layers, (1.0, 1.0)) == [8, 108]
  @test all(isnan, st_extract(layers, (4.0, 4.0)))
  @test st_extract(ra, points[1:2]; combine=vcat).value == [1, 5]

  # Same physical raster, with latitude stored south to north.
  ascending = rast(; A=reverse(ra.A; dims=2), b=ra.b,
                   lon=ra.lon, lat=reverse(ra.lat), cellsize=(1.0, 1.0))
  @test isequal(st_extract(ascending, points).value, result.value)

  # Do not snap nearby points onto a boundary.
  near = [(1 - 1e-8, 1.5), (1 + 1e-8, 1.5),
          (1.5, 1 - 1e-8), (1.5, 1 + 1e-8),
          (prevfloat(0.0), 1.5), (nextfloat(3.0), 1.5),
          (1.5, prevfloat(0.0)), (1.5, nextfloat(3.0))]
  @test isequal(vec(st_extract(ra, near).value), [4, 5, 8, 5, fill(NaN, 4)...])

  # Decimal extents: use bbox resolution, not differences of rounded centers.
  decimal = rast(reshape(1:9, 3, 3), bbox(100.1, 20.1, 100.4, 20.4))
  @test st_extract(decimal, [(100.2, 20.150000000000002),
                            (100.14999999999999, 20.199999999999999)]).value == [8 7]
  geographic = rast(reshape(1:48, 8, 6), bbox(110.7, 32.2, 110.72, 32.215))
  @test st_extract(geographic, (110.7025, 32.202500000000001)) == 41
end

@testset "extraction types and thin rasters" begin
  for T in (Float16, Float32, Float64, Int), bands in (1, 2)
    ra = rast(ones(T, 3, 3, bands), bbox(0, 0, 3, 3))
    @test eltype(st_extract(ra, [(0.5, 2.5)]).value) === T
    result = st_extract(ra, [(0.5, 2.5), (4.0, 4.0)])
    @test eltype(result.value) === float(T)
    @test all(==(one(T)), result.value[:, 1])
    @test all(isnan, result.value[:, 2])
    nodata = st_extract(ra, (4.0, 4.0))
    @test eltype(nodata) === float(T)
    @test all(isnan, nodata)
  end

  for dims in ((1, 1), (1, 3), (3, 1))
    ra = rast(reshape(1:prod(dims), dims), bbox(0, 0, 3, 3))
    points = [(x, y) for y in ra.lat for x in ra.lon]
    @test vec(st_extract(ra, points).value) == vec(ra.A)
    @test [st_extract(ra, p) for p in points] == vec(ra.A)
    @test st_extract(ra, (0.0, 3.0)) == first(ra.A)
    @test st_extract(ra, (3.0, 0.0)) == last(ra.A)
    @test isnan(st_extract(ra, (4.0, 4.0)))
  end
end
