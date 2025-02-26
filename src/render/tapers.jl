# Rendering of arbitrary segments
function to_polygons(f, len, s::Paths.TaperTrace; kwargs...)
    bnds = (zero(len), len)

    g = (t, sgn) -> begin
        d = Paths.direction(f, t) + sgn * 90.0°
        return f(t) + Paths.extent(s, t) * Point(cos(d), sin(d))
    end

    pgrid = adapted_grid(t -> Paths.direction(r -> g(r, 1), t), bnds; kwargs...)
    mgrid = adapted_grid(t -> Paths.direction(r -> g(r, -1), t), bnds; kwargs...)

    pts = [g.(mgrid, -1); @view (g.(pgrid, 1))[end:-1:1]]

    return Polygon(pts)
end

function to_polygons(f, len, s::Paths.TaperCPW; kwargs...)
    bnds = (zero(len), len)

    g =
        (t, sgn1, sgn2) -> begin
            d = Paths.direction(f, t) + sgn1 * 90.0°       # turn left (+) or right (-) of path
            offset = 0.5 * (Paths.gap(s, t) + Paths.trace(s, t))
            return f(t) + (sgn2 * 0.5 * Paths.gap(s, t) + offset) * Point(cos(d), sin(d))
        end

    ppgrid = adapted_grid(t -> Paths.direction(r -> g(r, 1, 1), t), bnds; kwargs...)
    pmgrid = adapted_grid(t -> Paths.direction(r -> g(r, 1, -1), t), bnds; kwargs...)
    mmgrid = adapted_grid(t -> Paths.direction(r -> g(r, -1, -1), t), bnds; kwargs...)
    mpgrid = adapted_grid(t -> Paths.direction(r -> g(r, -1, 1), t), bnds; kwargs...)

    ppts = [g.(pmgrid, 1, -1); @view (g.(ppgrid, 1, 1))[end:-1:1]]
    mpts = [g.(mpgrid, -1, 1); @view (g.(mmgrid, -1, -1))[end:-1:1]]

    return [Polygon(ppts), Polygon(mpts)]
end

# Optimized rendering of straight tapered segments
function to_polygons(segment::Paths.Straight{T}, s::Paths.TaperTrace; kwargs...) where {T}
    dir = direction(segment, zero(T))
    dp, dm = dir + 90.0°, dir - 90.0°

    ext_start = Paths.extent(s, zero(T))
    ext_end = Paths.extent(s, pathlength(segment))

    tangents = StaticArrays.@SVector [
        ext_start * Point(cos(dm), sin(dm)),
        ext_end * Point(cos(dm), sin(dm)),
        ext_end * Point(cos(dp), sin(dp)),
        ext_start * Point(cos(dp), sin(dp))
    ]

    a, b = segment(zero(T)), segment(pathlength(segment))
    origins = StaticArrays.@SVector [a, b, b, a]

    return Polygon(origins .+ tangents)
end

function to_polygons(segment::Paths.Straight{T}, s::Paths.TaperCPW; kwargs...) where {T}
    dir = direction(segment, zero(T))
    dp = dir + 90.0°

    ext_start = Paths.extent(s, zero(T))
    ext_end = Paths.extent(s, pathlength(segment))
    trace_start = Paths.trace(s, zero(T))
    trace_end = Paths.trace(s, pathlength(segment))

    tangents = StaticArrays.@SVector [
        Point(cos(dp), sin(dp)),
        Point(cos(dp), sin(dp)),
        Point(cos(dp), sin(dp)),
        Point(cos(dp), sin(dp))
    ]

    extents_p =
        StaticArrays.@SVector [0.5 * trace_start, 0.5 * trace_end, ext_end, ext_start]
    extents_m =
        StaticArrays.@SVector [ext_start, ext_end, 0.5 * trace_end, 0.5 * trace_start]

    a, b = segment(zero(T)), segment(pathlength(segment))
    origins = StaticArrays.@SVector [a, b, b, a]

    return [
        Polygon(origins .+ extents_p .* tangents),
        Polygon(origins .- extents_m .* tangents)
    ]
end
