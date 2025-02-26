function to_polygons(seg::Paths.Corner, ::Paths.SimpleTraceCorner; kwargs...)
    sgn = ifelse(seg.α >= 0.0°, 1, -1)

    ext = seg.extent * tan(sgn * seg.α / 2)
    p0 = seg.p0 - ext * Point(cos(seg.α0), sin(seg.α0))

    ∠A = seg.α0 + sgn * 90.0°
    p = Point(cos(∠A), sin(∠A))

    p1 = seg.extent * p + p0
    p2 = -seg.extent * p + p0
    p3 = p2 + 2ext * Point(cos(seg.α0), sin(seg.α0))
    p4 = p3 + 2ext * Point(cos(seg.α0 + seg.α), sin(seg.α0 + seg.α))

    arr = orientation(p1, p2, p3) > 0 ? [p1, p2, p3, p4] : [p4, p3, p2, p1]
    return Polygon(arr)
end
