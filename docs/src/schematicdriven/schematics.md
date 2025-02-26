# Schematics

## Schematic graphs

The schematic graph is a collection of logical information about the presence of components
and their connections. Its nodes are `ComponentNode`s, representing instances of `AbstractComponent`s
in the design. The edges are labeled with the names of the components' `Hook`s that are to
be fused together.

```@docs
SchematicDrivenLayout.SchematicGraph
SchematicDrivenLayout.ComponentNode
```

### Adding and connecting Components

```@docs
SchematicDrivenLayout.add_node!
SchematicDrivenLayout.fuse!
route!(::SchematicDrivenLayout.SchematicGraph,
    ::Paths.RouteRule,
    ::Pair{SchematicDrivenLayout.ComponentNode, Symbol},
    ::Pair{SchematicDrivenLayout.ComponentNode, Symbol},
    ::Any,
    ::Any)
SchematicDrivenLayout.RouteComponent
attach!(::SchematicDrivenLayout.SchematicGraph,
    ::S,
    ::Pair{T, Symbol},
    ::DeviceLayout.Coordinate
) where {S <: SchematicDrivenLayout.ComponentNode, T <: SchematicDrivenLayout.ComponentNode}
```

## Routing

A [`Paths.Route`](@ref) is like a `Path`, but defined implicitly in terms of its endpoints and rules (like "use only straight sections and 90 degree turns") for getting from one end to another. We can add `Route`s between components in a schematic using [`route!`](@ref), creating flexible connections that are only resolved after floorplanning has determined the positions of the components to be connected.

!!! info "Terminology note: connected component"
    
    There is a graph-theory term "connected component" unrelated to our `AbstractComponent`, indicating a subgraph whose nodes are all connected by edges and which isn't part of a larger such subgraph. For example, a fully connected graph has one connected component, the entire graph. Sometimes you may even hear these called "components", but below we use "connected component" for the graph-theory sense and "component" for `AbstractComponent`.

In more detail: Typically, each `fuse!` operation fully determines the position and rotation of the new node. For each connected component, the first node added to the `SchematicGraph` is positioned at the origin. A typical workflow could start by creating a node with a "chip template" component containing port [`hooks`](@ref SchematicDrivenLayout.hooks) and other boilerplate. We then `fuse!` any CPW launchers to the template. The devices in the middle of the chip are added and fused to one another without yet fixing their position relative to the chip template. Next, one connection is made between this connected component of devices and the connected component containing the template, fixing their relative positions. Since the chip template was added first, it will be centered at the origin.

At this point, further connections still need to be made between various device ports and CPW launchers positioned on the chip template, all of which are now fully constrained. Another constraint like those created by `fuse!` so far would overconstrain the layout, causing floorplanning with `plan` to fail. The remaining connections are instead made with [`route!`](@ref), which creates a `RouteNode` with edges to the `ComponentNode`s at its start and end points. Then, during `plan`, after the initial floorplanning phase has determined position of all fixed `Component`s, the route node finds a path between its start and end points.

!!! tip "Schematic connections without geometric constraints"
    
    You can add edges to the schematic graph that will be ignored during `plan` using the keyword `plan_skips_edge=true` in `fuse!`.

## Schematics

```@docs
SchematicDrivenLayout.Schematic
SchematicDrivenLayout.plan
SchematicDrivenLayout.check!
SchematicDrivenLayout.build!
SchematicDrivenLayout.render!(::SchematicDrivenLayout.AbstractCoordinateSystem, ::SchematicDrivenLayout.Schematic, ::SchematicDrivenLayout.LayoutTarget; kwargs...)
```

### Inspecting and manipulating schematics

Schematics allow you to produce layouts by specifying designs at a high level, in terms of
components and their connections. Often, this isn't quite enough to get the final artwork
you want. Maybe a wire needs to be routed differently, or a component needs to be modified
based on its position.

SchematicDrivenLayout allows you to inspect and modify a `Schematic` to make these sorts of
changes before `build` renders it into polygons.

When you construct a `SchematicGraph`, you don't need to know exactly where a component will
end up. You can usually calculate it yourself, but there are some built-in utilities to
simplify things.

```@docs
SchematicDrivenLayout.indexof(::SchematicDrivenLayout.ComponentNode, ::SchematicDrivenLayout.SchematicGraph)
SchematicDrivenLayout.bounds(::SchematicDrivenLayout.Schematic, ::SchematicDrivenLayout.ComponentNode)
SchematicDrivenLayout.center(::SchematicDrivenLayout.Schematic, ::SchematicDrivenLayout.ComponentNode)
SchematicDrivenLayout.hooks(::SchematicDrivenLayout.Schematic, ::SchematicDrivenLayout.ComponentNode)
SchematicDrivenLayout.origin(::SchematicDrivenLayout.Schematic, ::SchematicDrivenLayout.ComponentNode)
SchematicDrivenLayout.transformation(::SchematicDrivenLayout.Schematic, ::SchematicDrivenLayout.ComponentNode)
SchematicDrivenLayout.find_components
SchematicDrivenLayout.find_nodes
```

If you have a `ComponentNode`, you can replace its component using an arbitrary function
that takes the original component as an argument:

```@docs
SchematicDrivenLayout.replace_component!
```

You can also replace it (or all components of that type) using a function that takes the
position in schematic-global or wafer-global coordinates as an argument:

```@docs
SchematicDrivenLayout.position_dependent_replace!
```

One other useful trick allows a kind of interactive routing. When you view your final layout
built from the schematic, you may find that a route bends too sharply or goes too close to a
component. You can write down the points it needs to go to in the schematic's global coordinate system,
and add them as waypoints to the route. That is, if you go back to your layout script,
before you `build!` the layout, you can do something like

```julia
### original script
g = SchematicGraph("example")
...
route_node = route!(g, args...)
...
floorplan = plan(g)
check!(floorplan)
### modifications
# waypoints are global coordinates, not relative to the route's origin
route_node.component.global_waypoints = true
route_node.component.r.waypoints = [Point(600.0μm, -3000.0μm)]
route_node.component.r.waydirs = [90°]
### finalize
build!(floorplan)
```

Now the route in `route_node` is guaranteed to pass through the point (600.0μm, -3000.0μm)
on its way to its destination. If the `RouteRule`'s implementation uses `waydirs`, then it
will also have a direction of 90° at that point.

### Automatic crossover generation

You can automatically generate crossovers between `Path`s and `RouteComponent`s, including those nested within composite components. This is based on [`Path` intersection functionality](../paths.md#Intersections).

```@docs
SchematicDrivenLayout.crossovers!
```

### Checking schematics

It is often necessary to check that a planned `Schematic` obeys a set of constraints set
by the fabrication process. For instance, one may want to verify that all the junctions
in a floorplan are oriented in the right direction, e.g. pointing north. Instead of doing
this by eye, users should call `check!(sch::Schematic)`. In the future, this method could
run any number of checks, but at present it only checks the global orientation of
[_checkable_](@ref SchematicDrivenLayout.check_rotation) components via the following methods:

```@docs
SchematicDrivenLayout.rotations_valid
```

To be able to `build!` a floorplan (i.e. turn components into their `CoordinateSystem` geometries), users _must_
run `check!` first, otherwise `build` will throw an error.

## Visualization

SchematicDrivenLayout also contains some prototype schematic-visualization methods: if `GraphMakie` is present in the environment, then `GraphMakie.graphplot` can be used with a `SchematicGraph` to show nodes and edges or with `Schematic` to also place the components at their on-chip positions.
