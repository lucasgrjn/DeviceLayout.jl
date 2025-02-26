# Solid Models

A `Schematic` can be rendered to a `SolidModel`, creating a 3D geometry corresponding to the schematic, which can then be exported to a standard 3D format.
This is analogous to rendering the `Schematic` to a `Cell` in preparation for GDS export.

```@docs
render!(::DeviceLayout.SolidModel, ::SchematicDrivenLayout.Schematic, ::SchematicDrivenLayout.Target; kwargs...)
```

We use a `SolidModelTarget` to specify how a 2D geometry is rendered to 3D entities and physical groups.
This is analogous to the use of a `LayoutTarget` in `Cell` rendering, which specifies how to render entities with semantic metadata to `Cell` polygons and GDS layers.

```@docs
SchematicDrivenLayout.SolidModelTarget
```

Rendering using the `SolidModelTarget` applies the rendering option `solidmodel=true`. A pair of functions are provided for designating entities to be rendered or not based on that option (using `DeviceLayout.OptionalStyle`):

```@docs
SchematicDrivenLayout.not_solidmodel
SchematicDrivenLayout.only_solidmodel
```
