# Targets

```@docs
SchematicDrivenLayout.Target
SchematicDrivenLayout.LayoutTarget
SchematicDrivenLayout.ArtworkTarget
SchematicDrivenLayout.SimulationTarget
```

## Metadata handling

SchematicDrivenLayout provides the `facing` and `backing` functions to be used with a specific
interpretation of `level` in `SemanticMeta`. The level of a geometric entity describes the
vertical index of its substrate surface in a "flipchip"-style stack of substrates. Metadata types without a
level attribute will default to level 1.

```
▒    ...        ▒
▒   level 3 ↓   ▒
█████████████████
▒   level 2 ↑   ▒
▒               ▒
▒   level 1 ↓   ▒
█████████████████
▒   level 0 ↑   ▒
```

```@docs
SchematicDrivenLayout.backing
SchematicDrivenLayout.facing
```

## Rendering

A `Schematic` can be rendered to different geometry representations like `Cell` or `SolidModel` using different
`Target`s to control rendering options. See [`SchematicDrivenLayout.render!(::SchematicDrivenLayout.AbstractCoordinateSystem, ::SchematicDrivenLayout.Schematic, ::SchematicDrivenLayout.LayoutTarget)`](@ref) and [`SchematicDrivenLayout.render!(::DeviceLayout.SolidModel, ::SchematicDrivenLayout.Schematic, ::SchematicDrivenLayout.Target)`](@ref).

### Rendering flags

The built-in targets `ArtworkTarget` and `SimulationTarget` have the rendering options
`(artwork=true, simulation=false)` and `(artwork=false, simulation=true)`. A pair of
functions are provided for designating entities to be rendered or not based on the
`simulation` option (using `DeviceLayout.OptionalStyle`).

```@docs
SchematicDrivenLayout.not_simulated
SchematicDrivenLayout.only_simulated
```
