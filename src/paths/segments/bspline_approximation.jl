# Calculate approximation error
function _approximation_error(
    f::Paths.Segment{T},
    f_approx::Paths.BSpline,
    testvals=nothing
) where {T}
    # In this case, testvals are chosen based on f_approx, so they will be different every iteration
    tgrid, exact, exact_normal, _ = _testvals(f, f_approx)
    approx = f_approx.r.(tgrid)
    # f_approx may not be arclength-parameterized
    # but each exact point should be within 1nm of some line segment on
    # the discretization of f_approx
    idx_0 = 1
    maxerr = zero(T)
    for (p, normal) in zip(exact, exact_normal) # For each exact point, project it onto the approximate segments        
        # Find the first approx segment where the projection lies on that segment
        # (Starting with the approximate segment previously found)
        for idx = idx_0:(length(approx) - 1)
            intersects, ixn = Polygons.intersection(
                Polygons.Ray(p, p + normal),
                Polygons.LineSegment(approx[idx], approx[idx + 1])
            )
            if !intersects # Maybe we checked the ray in the wrong direction
                intersects, ixn = Polygons.intersection(
                    Polygons.Ray(p, p - normal),
                    Polygons.LineSegment(approx[idx], approx[idx + 1])
                )
            end
            if intersects # We found a projection that lies on the approximation
                idx_0 = idx
                maxerr = max(maxerr, norm(p - ixn))
                break
            end
        end
    end
    return maxerr
end

function _testvals(f::Paths.Segment{T}, f_approx::Paths.BSpline) where {T}
    l = Paths.pathlength(f)
    h(t) = Paths.Interpolations.hessian(f_approx.r, t)[1]
    tgrid = DeviceLayout.discretization_grid(h, _default_curve_atol(T))
    exact = f.(tgrid * l)
    dir = direction.(Ref(f), tgrid * l)
    normal = oneunit(T) * Point.(-sin.(dir), cos.(dir))
    return tgrid, exact, normal, nothing
end

function _testvals(
    f::Paths.ConstantOffset{T, BSpline{T}},
    f_approx::Paths.BSpline
) where {T}
    h(t) = Paths.Interpolations.hessian(f.seg.r, t)[1]
    tgrid = DeviceLayout.discretization_grid(h, _default_curve_atol(T))
    off = abs(getoffset(f))
    # Don't actually calculate the exact curve, we'll only check that the offset is correct
    exact = f.seg.r.(tgrid)
    tangent = first.(Paths.Interpolations.gradient.(Ref(f.seg.r), tgrid))
    normal = Point.(-gety.(tangent), getx.(tangent))
    return tgrid, exact, normal, off
end

function _testvals(f::Paths.GeneralOffset{T, BSpline{T}}, f_approx::Paths.BSpline) where {T}
    h(t) = Paths.Interpolations.hessian(f.seg.r, t)[1]
    tgrid = DeviceLayout.discretization_grid(h, _default_curve_atol(T))
    sgrid = t_to_arclength.(Ref(f.seg), tgrid)
    offsets = abs.(getoffset.(Ref(f), sgrid))
    # Don't actually calculate the exact curve, we'll only check that the offsets are correct
    exact = f.seg.r.(tgrid)
    tangent = first.(Paths.Interpolations.gradient.(Ref(f.seg.r), tgrid))
    normal = Point.(-gety.(tangent), getx.(tangent))
    return tgrid, exact, normal, offsets
end

# For BSpline offsets, use t parameterization rather than arclength (much faster)
function _approximation_error(
    f::Paths.ConstantOffset{T, BSpline{T}},
    f_approx::Paths.BSpline{T},
    testvals=_testvals(f, f_approx)
) where {T}
    tgrid, exact, exact_normal, off = testvals
    approx = f_approx.r.(tgrid)
    # f_approx may not be arclength-parameterized
    # but each exact point should be within 1nm of some line segment on
    # the discretization of f_approx
    idx_0 = 1
    maxerr = zero(T)
    for (p, normal) in zip(exact, exact_normal) # For each exact point, project it onto the approximate segments        
        # Find the first approx segment where the projection lies on that segment
        # (Starting with the approximate segment previously found)
        for idx = idx_0:(length(approx) - 1)
            intersects, ixn = Polygons.intersection(
                Polygons.Ray(p, p + normal),
                Polygons.LineSegment(approx[idx], approx[idx + 1])
            )
            if !intersects # Maybe we checked the ray in the wrong direction
                intersects, ixn = Polygons.intersection(
                    Polygons.Ray(p, p - normal),
                    Polygons.LineSegment(approx[idx], approx[idx + 1])
                )
            end
            if intersects # We found a projection that lies on the approximation
                idx_0 = idx
                maxerr = max(maxerr, abs(norm(p - ixn) - off))
                break
            end
        end
    end
    return maxerr
end

# For BSpline variable offsets, we need arclength-parameterized offsets for every test point
function _approximation_error(
    f::Paths.GeneralOffset{T, BSpline{T}},
    f_approx::Paths.BSpline{T},
    testvals=_testvals(f, f_approx)
) where {T}
    tgrid, exact, exact_normal, offsets = testvals
    approx = f_approx.r.(tgrid)
    # f_approx may not be arclength-parameterized
    # but each exact point should be within 1nm of some line segment on
    # the discretization of f_approx
    idx_0 = 1
    maxerr = zero(T)
    for (p, normal, off) in zip(exact, exact_normal, offsets) # For each exact point, project it onto the approximate segments        
        # Find the first approx segment where the projection lies on that segment
        # (Starting with the approximate segment previously found)
        for idx = idx_0:(length(approx) - 1)
            intersects, ixn = Polygons.intersection(
                Polygons.Ray(p, p + normal),
                Polygons.LineSegment(approx[idx], approx[idx + 1])
            )
            if !intersects # Maybe we checked the ray in the wrong direction
                intersects, ixn = Polygons.intersection(
                    Polygons.Ray(p, p - normal),
                    Polygons.LineSegment(approx[idx], approx[idx + 1])
                )
            end
            if intersects # We found a projection that lies on the approximation
                idx_0 = idx
                maxerr = max(maxerr, abs(norm(p - ixn) - off))
                break
            end
        end
    end
    return maxerr
end

function _sample_points(f::Paths.Segment{T}, num_points) where {T}
    s = range(zero(T), pathlength(f), length=num_points)
    return [p0(f), (f.(s[2:(end - 1)]))..., p1(f)]
end

function _sample_points(b::Paths.BSpline{T}, num_points) where {T}
    t = range(0.0, 1.0, length=num_points)
    return [p0(b), (b.r.(t[2:(end - 1)]))..., p1(b)]
end

function _sample_points(f::Paths.ConstantOffset{T, BSpline{T}}, num_points) where {T}
    b = f.seg
    t = range(0.0, 1.0, length=num_points)[2:(end - 1)]
    G = StaticArrays.@MVector [zero(Point{T})]
    perp(t) = begin
        Interpolations.gradient!(G, b.r, t)
        return Point(-G[1].y, G[1].x) / norm(G[1])
    end
    return [p0(f), (b.r.(t) .+ (f.offset * perp.(t)))..., p1(f)]
end

function _sample_points(f::Paths.GeneralOffset{T, BSpline{T}}, num_points) where {T}
    b = f.seg
    t = range(0.0, 1.0, length=num_points)[2:(end - 1)]
    G = StaticArrays.@MVector [zero(Point{T})]
    perp(t) = begin
        Interpolations.gradient!(G, b.r, t)
        return Point(-G[1].y, G[1].x) / norm(G[1])
    end
    return [
        p0(f),
        (b.r.(t) .+ (getoffset(f).(t_to_arclength.(Ref(b), t)) .* perp.(t)))...,
        p1(f)
    ]
end

function _initial_guess(f::Paths.Segment{T}) where {T}
    # Scaling tangents by pathlength/2 seems to improve convergence
    # (compared with scaling by pathlength)
    l = pathlength(f)
    t0 = Point(cos(α0(f)), sin(α0(f))) * l / 2
    t1 = Point(cos(α1(f)), sin(α1(f))) * l / 2
    return BSpline{T}([p0(f), p1(f)], t0, t1)
end

function _initial_guess(f::Paths.ConstantOffset{T, BSpline{T}}) where {T}
    # Original bspline has points evenly spaced in t from 0 to 1
    # If these were important for defining the original, the offset curve might need them
    b = f.seg
    pts = _sample_points(f, length(b.p))
    return BSpline{T}(pts, f.seg.t0, f.seg.t1) # tangents are the same
end

function _initial_guess(f::Paths.GeneralOffset{T, BSpline{T}}) where {T}
    # Original bspline has points evenly spaced in t from 0 to 1
    b = f.seg
    pts = _sample_points(f, length(b.p))
    # Non-constant offset means endpoint tangents may be different
    t0 = Point(cos(α0(f)), sin(α0(f))) * norm(f.seg.t0)
    t1 = Point(cos(α1(f)), sin(α1(f))) * norm(f.seg.t1)
    return BSpline{T}(pts, t0, t1)
end

_default_curve_atol(::Type{<:Real}) = 1e-3
_default_curve_atol(::Type{<:Length}) = 1nm

# Approximate f with a BSpline
function bspline_approximation(
    f::Paths.Segment{T};
    atol=_default_curve_atol(T),
    maxits=20
) where {T}
    # Sample points from f and use them to create the BSpline interpolation
    approx = _initial_guess(f)
    # Sample a dense set of points to test approximation against
    # (These testvals can be reused, although currently it's only reused for offsets of BSplines)
    testvals = _testvals(f, approx)
    # Calculate the maximum distance between approx and its projection onto the
    # discretization of the exact curve given by testvals
    err = _approximation_error(f, approx, testvals)
    # Double the number of interpolation points until the error is below tolerance
    refine = 1
    while err > atol
        refine > maxits && error(
            "Maximum iterations $maxits reached for B-spline approximation: Error $err > $atol"
        )
        # Double the number of interpolation segments with equal lengths in t
        npoints = 2 * length(approx.p) - 1
        approx.p = _sample_points(f, npoints)
        approx.t0 = approx.t0 / 2 # Halve the endpoint tangents when doubling the npoints
        approx.t1 = approx.t1 / 2
        Paths._update_interpolation!(approx)
        err = _approximation_error(f, approx, testvals)
        refine = refine + 1
    end
    return approx
end

bspline_approximation(b::Paths.BSpline; kwargs...) = copy(b)
