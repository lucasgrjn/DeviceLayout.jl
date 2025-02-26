## Styled geometry entities

Entities can also be "styled"—paired with a `GeometryEntityStyle`. This creates a `StyledEntity <: GeometryEntity` that still supports the [entity interface](geometry.md#Entities) (including the ability to be styled).

```@docs
DeviceLayout.GeometryEntityStyle
DeviceLayout.StyledEntity
DeviceLayout.entity
DeviceLayout.style(::DeviceLayout.StyledEntity)
DeviceLayout.styled
DeviceLayout.unstyled
DeviceLayout.unstyled_type
```

### Basic styles

The styles below can be applied to most entities for generic purposes.

```@docs
DeviceLayout.Plain
DeviceLayout.MeshSized
DeviceLayout.meshsized_entity
DeviceLayout.NoRender
DeviceLayout.OptionalStyle
DeviceLayout.optional_entity
DeviceLayout.ToTolerance
```

`AbstractPolygon`s also have the [`Polygons.Rounded`](@ref) style, while `ClippedPolygons` have [`Polygons.StyleDict`](@ref).
