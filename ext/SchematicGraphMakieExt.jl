module SchematicGraphMakieExt
using DeviceLayout, .SchematicDrivenLayout
import .SchematicDrivenLayout: indexof
import Graphs
import Graphs: edges, neighbors, src, dst, vertices
import MetaGraphs
import MetaGraphs:
    AbstractMetaGraph,
    MetaGraph,
    get_prop,
    has_prop,
    props,
    add_vertex!,
    add_edge!,
    rem_vertex!
import GraphMakie
import GraphMakie: graphplot
import Unitful: ustrip, nm

### Visualization with GraphPlot
"""
    graphplot(g::SchematicGraph)

Draw `g` using the `GraphMakie` package.

Nodes are labeled with the component name, while edges are labeled with the names
of the hooks on either endpoint (in arbitrary order).
"""
function GraphMakie.graphplot(g::SchematicGraph)
    nodelabel = ["[$(indexof(node, g))] $(node.id)" for node in nodes(g)]
    hook1s = get_prop.(g, edges(g), g[src.(edges(g))])
    hook2s = get_prop.(g, edges(g), g[dst.(edges(g))])
    edgelabel = ["$(string(h1)) -- $(string(h2))" for (h1, h2) in zip(hook1s, hook2s)]

    return graphplot(g, nlabels=nodelabel, elabels=edgelabel)
end

"""
    graphplot(sch::Schematic)

Draw the schematic graph of `sch`, attempting to arrange nodes according to the floorplan.

Nodes are labeled with the component name, while edges are labeled with the names
of the hooks on either endpoint (in arbitrary order).
"""
function GraphMakie.graphplot(sch::Schematic)
    g = sch.graph
    origins = (center.(Ref(sch), nodes(g)))
    origins = [round.(ustrip.(nm, p0)) for p0 in origins]
    xs = sort(unique(getx.(origins)))
    ys = sort(unique(gety.(origins)))

    spacing = 1.0
    grid = collect(Base.product(1:length(xs), 1:length(ys)))
    origins =
        [Point(grid[indexin(p0.x, xs), indexin(p0.y, ys)]...) * spacing for p0 in origins]

    # Space them out
    for v in vertices(g)
        if component(nodes(g)[v]) isa RouteComponent
            origins[v] =
                sum([origins[v2] for v2 in neighbors(g, v)]) / length(neighbors(g, v))
        else
            for v2 in neighbors(g, v)
                if (sum((origins[v] - origins[v2]) .^ 2)) < 0.4
                    dir = in_direction(hooks(sch, g[v2], get_prop(g, g[v], g[v2], g[v2])))
                    origins[v2] = origins[v2] + 0.5 * Point(cos(dir), sin(dir))
                end
            end
        end
    end

    nodelabel = ["[$(indexof(node, g))] $(node.id)" for node in nodes(g)]
    hook1s = get_prop.(g, edges(g), g[src.(edges(g))])
    hook2s = get_prop.(g, edges(g), g[dst.(edges(g))])
    edgelabel = ["$(string(h1)):$(string(h2))" for (h1, h2) in zip(hook1s, hook2s)]

    return graphplot(g, layout=(_) -> origins, nlabels=nodelabel, elabels=edgelabel)
end

MIMETypes = Union{
    MIME"image/png",
    MIME"image/svg+xml",
    MIME"application/pdf",
    MIME"application/postscript"
}
function Base.show(io, mime::MIMETypes, g::SchematicGraph; options...)
    return show(io, mime, graphplot(g); options...)
end
function Base.show(io, mime::MIMETypes, sch::Schematic; options...)
    return show(io, mime, graphplot(sch); options...)
end

end
