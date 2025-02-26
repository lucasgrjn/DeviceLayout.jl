# Functions useful for ExamplePDK but not suitable/ready for base DeviceLayout/SchematicDrivenLayout

"""
    filter_params(subcomp, comp::AbstractComponent; except=Symbol[])

Return the parameters of `comp` that share a name with parameters of `subcomp`.

Intended as a utility for passing parameters down to subcomponents.
`subcomp` can be an `AbstractComponent` instance or subtype.
"""
function filter_params(subcomp, comp; except=Symbol[])
    return filter(
        kv -> first(kv) in parameter_names(subcomp) && !(first(kv) in except),
        pairs(parameters(comp))
    )
end

"""
    tap!(path::Path, sty::Paths.SimpleCPW=laststyle(path); location=1)

Generate a new path branching off from an initial path.

Location should be `1` for a right-hand tap and `-1` for a left-hand tap.

To illustrate, we start with this input `path`, where the double arrow indicates
the forward direction from the endpoint:

    input path
        ⇑
    ███   ███
    ███   ███
    ███   ███

Then after `tap!`, the input path is modified and an output
path is returned as follows (with `location == 1`):

    input path
        ⇑
    ███   ███ ↕ sty.gap
    ███                            ⤒
    ███      ⇒ output path     sty.trace
    ███                            ⤓
    ███   ███ ↕ sty.gap
    ███   ███
    ███   ███
    ███   ███
"""
function tap!(path::Path, sty::Paths.SimpleCPW=laststyle(path); location=1)
    main_sty = laststyle(path)
    main_sty isa Paths.CPW || error("Last path style should be CPW for a CPW tap")
    main_trace = Paths.trace(main_sty, pathlength(path[end]))
    main_gap = Paths.gap(main_sty, pathlength(path[end]))
    # Add a virtual segment to get (next to) tap start
    norender = Paths.SimpleNoRender(main_trace + main_gap, virtual=true)
    straight!(path, Paths.extent(sty), norender)
    # Create tap path
    tap = Path(
        p1(path) - sign(location) * main_trace / 2 * Point(-sin(α1(path)), cos(α1(path))),
        α0=α1(path) - sign(location) * 90°,
        metadata=path.metadata,
        name=uniquename("tap")
    )
    straight!(tap, main_gap, sty)
    # Extend virtual segment the rest of the way
    straight!(path, Paths.extent(sty), norender)
    # Attach a rectangle to fill in the gap opposite the tap
    cs = CoordinateSystem(uniquename("tap_cut"), nm)
    place!(
        cs,
        centered(Rectangle(2 * Paths.extent(sty), laststyle(path).gap)),
        path.metadata
    )
    attach!(path, sref(cs), zero(main_gap), location=(-location))

    return tap
end

"""
    bridge_geometry(style::Paths.SimpleCPW)

Return a `CoordinateSystem` with a simple scaffolded bridge that spans `style`.
"""
function bridge_geometry(style::Paths.SimpleCPW)
    cs = CoordinateSystem(uniquename("bridge"))
    h_ground_ground = 2 * Paths.extent(style)
    place!(cs, centered(Rectangle(10μm, h_ground_ground + 20μm)), LayerVocabulary.BRIDGE)
    place!(
        cs,
        centered(Rectangle(16μm, h_ground_ground + 10μm)),
        LayerVocabulary.BRIDGE_BASE
    )
    return cs
end

"""
    add_bridges!(schematic, bridge=FEEDLINE_BRIDGE; spacing=200μm, margin=50μm)

Example utility for adding bridges. Not optimized for microwave properties.

Finds all top-level `Path`s and `RouteComponent`s in `schematic`.
For `Path`s, places a bridge in the middle of any path segment with length of at least `margin`.
For `RouteComponent`s, places a bridge at every `spacing`, with no bridges within a margin of
`margin` from the start and end.
"""
function add_bridges!(schematic, bridge; spacing=500μm, margin=50μm)
    ref = sref(bridge)
    for idx in find_components(Path, schematic.graph, depth=1)
        path = component(schematic.graph[idx])
        contains(path.name, "launcher") && continue
        add_bridges!(path, bridge; margin)
    end
    for idx in find_components(RouteComponent, schematic.graph, depth=1)
        routecomp = component(schematic.graph[idx])
        path = SchematicDrivenLayout.path(routecomp)
        attach!(routecomp, ref, margin:spacing:(pathlength(path) - margin))
    end
end

"""
    add_bridges!(path::Path, bridge; margin=50μm)

Example utility for adding bridges. Not optimized for microwave properties.

Places a bridge in the middle of any path segment with length of at least `margin`.
"""
function add_bridges!(path::Path, bridge; margin=50μm)
    isnothing(bridge) && return
    ref = sref(bridge)
    for (i, pathnode) in enumerate(path)
        Paths.undecorated(pathnode.sty) isa Paths.SimpleCPW || continue
        pathlength(pathnode) < margin && continue
        attach!(path, ref, pathlength(pathnode) / 2, i=i)
    end
end
