"""
    abstract type OffsetSegment{T, S <: Segment{T}} <: ContinuousSegment{T}

Supertype for "offset curves", defined by a segment with a signed offset.

An offset curve follows the original curve with some offset in the normal direction,
which may be constant (`ConstantOffset`) or vary with pathlength (`GeneralOffset`).
The curve normal is rotated +90° from the curve tangent, so that a positive offset
corresponds to the "left-hand side" of the curve.

Note that `pathlength(seg::OffsetSegment)` is the length of the original curve, and that
the offset segment is parameterized by the arclength of the original curve.

Used internally during SolidModel rendering. Not appropriate for use in Paths.
"""
abstract type OffsetSegment{T, S <: Segment{T}} <: ContinuousSegment{T} end

"""
    struct ConstantOffset{T,S} <: OffsetSegment{T,S}
"""
struct ConstantOffset{T, S} <: OffsetSegment{T, S}
    seg::S
    offset::T
end

"""
    struct ConstantOffset{T,S} <: OffsetSegment{T,S}
"""
struct GeneralOffset{T, S, U} <: OffsetSegment{T, S}
    seg::S
    offset::U
end

copy(s::OffsetSegment) = OffsetSegment(copy(s.seg), s.offset)
getoffset(s::ConstantOffset, l...) = s.offset
getoffset(s::GeneralOffset{T}) where {T} = l -> uconvert(unit(T), s.offset(l))
getoffset(s::GeneralOffset, l) = getoffset(s)(l)

# "Pathlength" is length of underlying segment!
# That's one reason OffsetSegments are not part of the public API.
pathlength(s::OffsetSegment) = pathlength(s.seg)
# You could get the "true" length using something like this:
# Assuming that s.seg has less than 1 winding!
# pathlength(s.seg) - s.offset * uconvert(NoUnits, α1(s.seg) - α0(s.seg))
# (See _curvelength)
# But this way offset segment is simply parameterized by arclength of the underlying curve
(s::ConstantOffset)(t) =
    s.seg(t) + s.offset * Point(-sin(direction(s.seg, t)), cos(direction(s.seg, t)))
(s::GeneralOffset{T})(t) where {T} =
    s.seg(t) + getoffset(s, t) * Point(-sin(direction(s.seg, t)), cos(direction(s.seg, t)))

p0(s::ConstantOffset) = p0(s.seg) + s.offset * Point(-sin(α0(s.seg)), cos(α0(s.seg)))
p1(s::ConstantOffset) = p1(s.seg) + s.offset * Point(-sin(α1(s.seg)), cos(α1(s.seg)))
α0(s::ConstantOffset) = α0(s.seg)
α1(s::ConstantOffset) = α1(s.seg)
direction(s::ConstantOffset, t) = direction(s.seg, t)

p0(s::GeneralOffset{T}) where {T} =
    p0(s.seg) + getoffset(s, zero(T)) * Point(-sin(α0(s.seg)), cos(α0(s.seg)))
α0(s::GeneralOffset{T}) where {T} = direction(s, zero(T))
direction(s::GeneralOffset, t) =
    direction(s.seg, t) + atan(ForwardDiff.derivative(s.offset, t))

summary(s::OffsetSegment) = summary(s.seg) * " offset by $(s.offset)"

# Methods for creating offset segments
OffsetSegment(seg::S, offset::Coordinate) where {T, S <: Segment{T}} =
    ConstantOffset{T, S}(seg, offset)
OffsetSegment(seg::S, offset::U) where {T, S <: Segment{T}, U} =
    GeneralOffset{T, S, U}(seg, offset)
offset(seg::Segment, s) = OffsetSegment(seg, s)
offset(seg::ConstantOffset, s::Coordinate) = offset(seg.seg, s + seg.offset)
offset(seg::ConstantOffset, s) = offset(seg.seg, t -> s(t) + seg.offset)
offset(seg::GeneralOffset, s::Coordinate) = offset(seg.seg, t -> s + seg.offset(t))
offset(seg::GeneralOffset, s) = offset(seg.seg, t -> s(t) + seg.offset(t))

function transform(x::ConstantOffset, f::Transformation)
    y = deepcopy(x)
    xrefl(f) && change_handedness!(y)
    setα0p0!(y.seg, rotated_direction(α0(y.seg), f), f(p0(y.seg)))
    return y
end

# Define outer constructors for Turn and Straight from
Straight(x::ConstantOffset{T, Straight{T}}) where {T} =
    Straight{T}(x.seg.l, p0=p0(x), α0=α0(x))
function Turn(x::ConstantOffset{T, Turn{T}}) where {T}
    return Turn(
        x.seg.α,
        x.seg.r + (abs(x.seg.α) > x.seg.α ? x.offset : -x.offset),
        p0=p0(x),
        α0=x.seg.α0
    )
end
# Methods for true length of offset curves
# Note that t parameterization is not necessarily arclength parameterization for BSplines
# (Or for offsets of BSplines)
function check_interval(f, t0, t1)
    valid = (t0 >= 0 && t1 <= 1 && t0 < t1)
    return valid || error(
        "Invalid interval [$t0, $t1] on $f. Interval should be of the form [t0, t1] with 0 <= t0 < t1 <= 1."
    )
end

function arclength(f::Segment, t1::Real=1.0; t0::Real=0.0)
    check_interval(f, t0, t1)
    return pathlength(f) * (min(t1, 1.0) - max(t0, 0.0))
end

function arclength(f::OffsetSegment{T}, t1::Real=1.0; t0::Real=0.0) where {T}
    check_interval(f, t0, t1)
    l_tot = pathlength(f)
    return quadgk(t -> norm(Paths.ForwardDiff.derivative(f, t)), t0 * l_tot, t1 * l_tot)[1]
end

# For underlying BSplines we can calculate the excess length due to curvature directly
function arclength(f::ConstantOffset{T, BSpline{T}}, t1::Real=1.0; t0::Real=0.0) where {T}
    check_interval(f, t0, t1)
    l_orig = arclength(f.seg, t1; t0=t0)
    G = StaticArrays.@MVector [zero(Point{T})]
    H = StaticArrays.@MVector [zero(Point{T})]
    l_extra = f.offset * quadgk(t -> _curvature_arclength!(f.seg, G, H, t), t0, t1)[1]
    return l_orig + l_extra
end

# For general offsets we also have to add the contribution from varying offset
function arclength(f::GeneralOffset{T, BSpline{T}}, t1::Real=1.0; t0::Real=0.0) where {T}
    check_interval(f, t0, t1)
    # If we unwind the curvature, each offset segment will have length (sqrt(1+x^2) ds)
    # where x is the derivative of the offset
    l_orig = quadgk(
        s -> sqrt(1 + (ForwardDiff.derivative(f.offset, s))^2),
        t_to_arclength(f.seg, t0),
        t_to_arclength(f.seg, t1)
    )[1]
    # Now add curvature contribution as usual
    G = StaticArrays.@MVector [zero(Point{T})]
    H = StaticArrays.@MVector [zero(Point{T})]
    l_extra = quadgk(
        t ->
            _curvature_arclength!(f.seg, G, H, t) * getoffset(f, t_to_arclength(f.seg, t)),
        t0,
        t1
    )[1]
    return l_orig + l_extra
end

function _curvature_arclength!(b::BSpline, G, H, t)
    # Curvature for arbitrary parameterization divides by another norm(g)
    # But when integrating over arclength we need to multiply by norm(g)
    Paths.Interpolations.gradient!(G, b.r, t)
    Paths.Interpolations.hessian!(H, b.r, t)
    return uconvert(NoUnits, (G[1].y * H[1].x - G[1].x * H[1].y) / norm(G[1])^2)
end
