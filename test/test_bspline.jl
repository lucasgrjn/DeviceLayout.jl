@testset "BSpline" begin
    po0 = Point(1.0μm, 1.0μm)
    po1 = Point(1000.0μm, -20.0μm)

    g0 = 100 * Point(1.0μm, 1.0μm) / sqrt(2)
    g1 = 100 * Point(0.6μm, -0.8μm)

    b = Paths.BSpline([po0, po1], g0, g1)

    # Segment properties
    @test Paths.arclength_to_t(b, Paths.t_to_arclength(b, 0.6)) ≈ 0.6 rtol = 1e-9
    @test Paths.Interpolations.gradient(b.r, 0.0)[1] ≈ g0 rtol = 1e-9
    @test Paths.Interpolations.gradient(b.r, 1.0)[1] ≈ g1 rtol = 1e-9
    @test Paths.curvatureradius(b, Paths.t_to_arclength(b, 0.0)) < zero(1.0μm)
    @test Paths.curvatureradius(b, Paths.t_to_arclength(b, 1.0)) < zero(1.0μm)

    # Reflect about y = 0, check the curvature radius changes
    g0 = 100 * Point(1.0μm, -1.0μm) / sqrt(2)
    g1 = 100 * Point(0.6μm, 0.8μm)
    b2 = Paths.BSpline([Point(1.0μm, -0.0μm), Point(1000.0μm, 20.0μm)], g0, g1)

    @test Paths.arclength_to_t(b, Paths.t_to_arclength(b, 0.6)) ≈ 0.6 rtol = 1e-9
    @test Paths.Interpolations.gradient(b2.r, 0.0)[1] ≈ g0 rtol = 1e-9
    @test Paths.Interpolations.gradient(b2.r, 1.0)[1] ≈ g1 rtol = 1e-9
    @test Paths.curvatureradius(b2, Paths.t_to_arclength(b2, 0.0)) > zero(1.0μm)
    @test Paths.curvatureradius(b2, Paths.t_to_arclength(b2, 1.0)) > zero(1.0μm)

    # Extending a Path
    path1 = Path(po0, α0=90°)
    bspline!(
        path1,
        [Point(1.0μm, 1001.0μm)],
        90°,
        Paths.SimpleCPW(20μm, 10μm),
        endpoints_speed=1.0μm
    )
    @test pathlength(path1) ≈ 1000μm rtol = 1e-9
    @test Paths.Interpolations.gradient(path1.nodes[1].seg.r, 0.0)[1] ≈ Point(0.0μm, 1.0μm) rtol =
        1e-9
    @test p1(path1) ≈ Point(1.0μm, 1001.0μm) rtol = 1e-9

    bspline!(
        path1,
        [Point(200.0μm, 500.0μm), po1],
        45°,
        Paths.TaperCPW(20μm, 10μm, 10μm, 5μm),
        endpoints_speed=100.0μm
    )
    @test Paths.Interpolations.gradient(path1.nodes[2].seg.r, 0.0)[1] ≈
          Point(0.0μm, 100.0μm) rtol = 1e-9
    @test Paths.Interpolations.gradient(path1.nodes[2].seg.r, 1.0)[1] ≈
          Point(100.0μm, 100.0μm) / sqrt(2) rtol = 1e-9

    # Constructor
    @test_throws ErrorException Paths.BSpline{Int}([Point(1, 1)], Point(2, 2), Point(3, 3))
    b2 = Paths.BSpline(
        [Point(1, 1), Point(1000, -20)],
        Point(1 / sqrt(2), 1 / sqrt(2)),
        Point(6, -8)
    )
    @test eltype(b2) == Float64

    # Convert
    b3 = convert(Paths.BSpline{typeof(1.0mm)}, b)
    @test eltype(b3) == typeof(1.0mm)

    # Split
    a1, a2 = Paths._split(b, 500μm)
    tsplit = Paths.arclength_to_t(b, 500μm)
    @test p0(a1) == po0
    @test p1(a1) == p0(a2)
    @test p1(a2) == po1
    @test a1(Paths.t_to_arclength(a1, 0.7)) ≈ b.r(tsplit * 0.7) atol = 1e-6 * μm
    @test a2(Paths.t_to_arclength(a2, 0.2)) ≈ b.r(tsplit + (1 - tsplit) * 0.2) atol =
        1e-9 * μm
    @test pathlength(a1) ≈ 500μm rtol = 1e-9
    @test pathlength(a2) ≈ pathlength(b) - 500μm rtol = 1e-9

    # Splice
    b4 = Paths.BSpline(
        [Point(-100, -200), Point(200, 200), Point(100, -300), Point(500, 0)],
        Point(-100, 800),
        Point(500, 100)
    )

    pa2 = Paths.split(Paths.Node(b4, Paths.Trace(10)), 500)
    tsplit2 = Paths.arclength_to_t(b4, 500)
    a4 = pa2[end].seg
    c = Cell{Float64}("bsp")
    render!(c, pa2, GDSMeta(); atol=0.1)

    # Prepare manual splice transform for comparison
    translate = Translation(p0(b4) - p0(a4))
    rotate = Rotation(α0(b4) - α0(a4))
    splice_transform = Translation(p0(b4)) ∘ rotate ∘ Translation(-p0(b4)) ∘ translate

    Paths.splice!(pa2, 1)
    @test p0(pa2) == p0(b4)
    @test α0(pa2) == α0(b4)
    # Check an arbitrary point to make sure we have just rotated and translated a curve segment
    @test (pa2[1].seg)(Paths.t_to_arclength(pa2[1].seg, 0.2)) ≈
          splice_transform(b4.r(tsplit2 + (1 - tsplit2) * 0.2)) atol = 1e-9
end
