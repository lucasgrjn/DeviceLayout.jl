```@meta
DocTestSetup = quote
    using Unitful, DeviceLayout
    using Unitful: °
end
```

# Coordinate Transformations

The mechanism for affine transformations is largely provided by the
[`CoordinateTransformations.jl`](https://github.com/FugroRoames/CoordinateTransformations.jl)
package. For convenience, the documentation for `Translation` and `compose` is
reproduced below from that package. We also provide convenience constructors for `Reflection`
across a specified line and `Rotation` around a specified point, as well as a `ScaledIsometry`
type that represents transformations restricted to those that preserve angles.

## Creating transformations

```@docs
    CoordinateTransformations.compose
    CoordinateTransformations.Translation
    Reflection
    XReflection
    YReflection
    Rotation
    RotationPi
    ScaledIsometry
```

Transformations can also be inverted with `inv`.

## Applying transformations

Coordinate transformations can be applied to any `AbstractGeometry` object, creating a new
object with its coordinates transformed. Transformations created using the constructors in
the above section can be applied directly to objects like a function. Here's an example with a `Rectangle`:

```jldoctest
julia> r = Rectangle(2, 1)
Rectangle{Int64}((0,0), (2,1))

julia> trans = Translation(10, 10)
Translation(10, 10)

julia> trans = Rotation(90°) ∘ trans
AffineMap([0.0 -1.0; 1.0 0.0], [-10.0, 10.0])

julia> trans(r)
Rectangle{Float64}((-11.0,10.0), (-10.0,12.0))
```

### Simple transformations

There are methods for conveniently applying simple transformations:

```@docs
    centered
    magnify
    reflect_across_line
    reflect_across_xaxis
    rotate
    rotate90
    translate
    +(::DeviceLayout.AbstractGeometry, ::Point)
    -(::DeviceLayout.AbstractGeometry, ::Point)
    *(::DeviceLayout.AbstractGeometry, a::Real)
    /(::DeviceLayout.AbstractGeometry, a::Real)
```

### Alignment

There are also methods to apply transformations that align objects using the edges of their
bounding boxes.

```@docs
    Align.above
    Align.below
    Align.leftof
    Align.rightof
    Align.flushbottom
    Align.flushtop
    Align.flushleft
    Align.flushright
    Align.centered_on
    Align.aligned_to
```

## Inspecting transformations

```@docs
    isapprox_angle
    isapprox_cardinal
    mag
    origin
    preserves_angles
    rotated_direction
    rotation
    DeviceLayout.Transformations.rounding_safe
    xrefl
```

## Implementation details

Geometry objects implement specializations of the `transform` function that determine how
they behave under transformations:

```@docs
    transform
```

This allows special handling for certain types paired with certain transformations.
For example, a `Rectangle` is by definition axis-aligned. If a transformed `Rectangle`
would still be axis-aligned (for example, the result of a translation or 90° rotation),
the result will still be a `Rectangle`, as in the example above; otherwise, it will be a `Polygon`.

"Chiral" geometry objects (those that can be either left- or right-handed)
can also implement special handling for transformations that include a reflection (which
changes handedness).
