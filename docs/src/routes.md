# Routes

A [`Paths.Route`](@ref) implicitly defines a `Path` between two points, following a
[`Paths.RouteRule`](@ref), with a specified start and end direction. It's a building block
for "interactive autorouting". For simple cases, it lets you define a path between two
points without having to do any geometric calculations yourself. In more complicated cases,
you can provide additional waypoint constraints to guide it.

To draw a `Path` from a `Route`, you can use [`Path(r::Route, sty)`](@ref).

More often, you won't work directly with the `Route` type but will instead use [`route!`](@ref) to extend an existing path to an endpoint according to the rules you specify.

## Route API

```@docs
    Paths.Route
```

### Route rules

```@docs
    Paths.RouteRule
    Paths.BSplineRouting
    Paths.StraightAnd90
    Paths.StraightAnd45
    Paths.CompoundRouteRule
```

### Route drawing

Calling `Path(r::Route, sty::Style)` creates a new `Path` at `p0(r)`, then extends it to
`p1(r)` using `route!`. The default implementation of `route!`, for a generic `RouteRule`,
first calls `reconcile!` to validate and modify waypoints as necessary. It then calls
`_route!` to draw the path, which by default calls `_route_leg!` to each waypoint in
order. (A "leg" is an informal abstraction describing the "unit" that the `RouteRule`
works with, which may be more than one `Paths.Segment`. For example, each leg in
`StraightAnd90` routing has a single 90-degree bend with up to one straight segment on
either side, or a single straight segment if no bend is necessary.)

```@docs
    Paths.Path(::Paths.Route, ::Paths.Style)
    Paths.route!
    Paths.reconcile!(::Paths.Path, ::Point, ::Any, ::Paths.RouteRule, ::Any, ::Any)
```

### Route inspection

A `Route` supports endpoint inspection much like a `Path` does:

```@docs
    Paths.p0(::Paths.Route)
    Paths.α0(::Paths.Route)
    Paths.p1(::Paths.Route)
    Paths.α1(::Paths.Route)
```
