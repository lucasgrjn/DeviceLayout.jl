function to_polygons(f, len, s::Paths.Strands; kwargs...)
    bnds = (zero(len), len)

    p = Polygon{typeof(len)}[]
    g =
        (t, sgn1, idx, sgn2) -> begin
            d = Paths.direction(f, t) + sgn1 * 90.0°       # turn left (+) or right (-) of path
            offset = Paths.offset(s, t) + Paths.width(s, t) / 2
            strand_offset = idx * (Paths.spacing(s, t) + Paths.width(s, t))
            return f(t) +
                   (sgn2 * Paths.width(s, t) / 2 + offset + strand_offset) *
                   Point(cos(d), sin(d))
        end
    for i = 0:(Paths.num(s) - 1)
        ppgrid = adapted_grid(t -> Paths.direction(r -> g(r, 1, i, 1), t), bnds; kwargs...)
        pmgrid = adapted_grid(t -> Paths.direction(r -> g(r, 1, i, -1), t), bnds; kwargs...)
        mmgrid =
            adapted_grid(t -> Paths.direction(r -> g(r, -1, i, -1), t), bnds; kwargs...)
        mpgrid = adapted_grid(t -> Paths.direction(r -> g(r, -1, i, 1), t), bnds; kwargs...)

        ppts = [g.(pmgrid, 1, i, -1); @view (g.(ppgrid, 1, i, 1))[end:-1:1]]
        mpts = [g.(mpgrid, -1, i, 1); @view (g.(mmgrid, -1, i, -1))[end:-1:1]]

        push!(p, Polygon(uniquepoints(ppts)))
        push!(p, Polygon(uniquepoints(mpts)))
    end
    return p
end

function to_polygons(
    segment::Paths.Straight{T},
    s::Paths.SimpleStrands;
    kwargs...
) where {T}
    dir = direction(segment, zero(T))
    dp = dir + 90.0°

    tangents = [
        Point(cos(dp), sin(dp)),
        Point(cos(dp), sin(dp)),
        Point(cos(dp), sin(dp)),
        Point(cos(dp), sin(dp))
    ]

    p = Polygon{T}[]
    for i = 0:(Paths.num(s) - 1)
        i_offset = i * (Paths.width(s) + Paths.spacing(s))
        o = Paths.offset(s) + i_offset
        ow = o + Paths.width(s)

        ext = Paths.extent(s, zero(T))

        extents_p = [o, o, ow, ow]
        extents_m = [ow, ow, o, o]

        a, b = segment(zero(T)), segment(pathlength(segment))
        origins = [a, b, b, a]

        push!(p, Polygon(origins .+ extents_p .* tangents))
        push!(p, Polygon(origins .- extents_m .* tangents))
    end
    return p
end
