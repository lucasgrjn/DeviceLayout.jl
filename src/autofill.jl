module Autofill

# Autofill, halos, and related functions
using ..Cells
using ..CoordinateSystems
using ..Polygons
using ..Paths
import ..Polygons: AbstractPolygon
import ..CoordinateSystems: AbstractCoordinateSystem
import DeviceLayout: Coordinate, GDSMeta, GeometryStructure, Meta, Point, SemanticMeta
import DeviceLayout:
    element_metadata,
    layer,
    layer_inclusion,
    place!,
    render!,
    halo,
    flatten,
    transformation,
    sref

export autofill!, make_halo

"""
    make_halo(delta, inner_delta=nothing; only_layers=[], ignore_layers=[])

Returns a function `(c)->halo(c, delta, inner_delta; only_layers, ignore_layers)`
for generating `GeometryStructure` halos.

Any entities in layers in `ignore_layers` will be skipped.
If `only_layers` is not empty, only those layers will be used to generate the halo.
Layers for inclusion and exclusion can be provided as layer name `Symbol`s, in which case
only the layer name needs to be matched, or as full `DeviceLayout.Meta` objects, in which case all
metadata fields (e.g., index and level for `SemanticMeta`) must match.
"""
function make_halo(delta, inner_delta=nothing; only_layers=[], ignore_layers=[])
    return (c) ->
        halo(c, delta, inner_delta; only_layers=only_layers, ignore_layers=ignore_layers)
end

"""
    autofill!(cs::AbstractCoordinateSystem,
        filler_cs::AbstractCoordinateSystem,
        grid_x::AbstractArray,
        grid_y::AbstractArray,
        exclusion)

Add references to `filler_cs` inside `cs` at grid points not in any `exclusion` polygon.

The `exclusion` argument may be

  - a `Coordinate` denoting an offset used to generate the excluded region from shapes in `cs`
  - a `CoordinateSystem` or `Cell` containing the geometry of the excluded region
  - a `Function` creating a `CoordinateSystem` or `Cell` from `cs`
  - an `AbstractArray{<:AbstractPolygon}`

Returns the origins of references.
"""
function autofill!(
    cs::AbstractCoordinateSystem,
    filler_cs::AbstractCoordinateSystem,
    grid_x::AbstractArray,
    grid_y::AbstractArray,
    exclusion::AbstractArray{<:AbstractPolygon}
)
    in_poly = gridpoints_in_polygon(exclusion, grid_x, grid_y)

    outpoly = findall((!), in_poly)
    origins = map(outpoly) do idx
        return Point(grid_x[first(idx.I)], grid_y[last(idx.I)])
    end

    addref!.(cs, filler_cs, origins)
    return origins
end

autofill!(
    cs::AbstractCoordinateSystem,
    filler_cs::AbstractCoordinateSystem,
    grid_x::AbstractArray,
    grid_y::AbstractArray,
    exclusion_delta::Coordinate
) = autofill!(cs, filler_cs, grid_x, grid_y, halo(cs, exclusion_delta))

autofill!(
    cs::AbstractCoordinateSystem,
    filler_cs::AbstractCoordinateSystem,
    grid_x::AbstractArray,
    grid_y::AbstractArray,
    exclusion_method::Function
) = autofill!(cs, filler_cs, grid_x, grid_y, exclusion_method(cs))

autofill!(
    cs::AbstractCoordinateSystem,
    filler_cs::AbstractCoordinateSystem,
    grid_x::AbstractArray,
    grid_y::AbstractArray,
    exclusion_cs::AbstractCoordinateSystem
) = autofill!(cs, filler_cs, grid_x, grid_y, Cell(exclusion_cs, map_meta=(_) -> GDSMeta()))

autofill!(
    cs::AbstractCoordinateSystem,
    filler_cs::AbstractCoordinateSystem,
    grid_x::AbstractArray,
    grid_y::AbstractArray,
    exclusion_cell::Cell
) = autofill!(cs, filler_cs, grid_x, grid_y, elements(flatten(exclusion_cell)))

"""
    halo(cs::CoordinateSystem, outer_delta, inner_delta=nothing; only_layers=[],
        ignore_layers=[])
    halo(cs::Cell, outer_delta, inner_delta=nothing; only_layers=[], ignore_layers=[])

A coordinate system of type `typeof(cs)` with halos for all entities in `cs`, tracing through `cs.refs`.

Any entities in layers in `ignore_layers` will be skipped.
If `only_layers` is not empty, only those layers will be used to generate the halo.
Layers for inclusion and exclusion can be provided as layer name `Symbol`s, in which case
only the layer name needs to be matched, or as full `DeviceLayout.Meta` objects, in which case all
metadata fields (e.g., index and level for `SemanticMeta`) must match.

The orientations of polygons must be consistent, such that outer polygons share the same
orientation, and any holes have the opposite orientation. Additionally, any holes should be
contained within outer polygons; offsetting hole edges may create positive artifacts at
corners.
"""
function halo(
    cs::CoordinateSystem{S},
    outer_delta,
    inner_delta=nothing;
    only_layers=[],
    ignore_layers=[],
    memoized_halos=Dict{GeometryStructure, GeometryStructure}()
) where {S}
    haskey(memoized_halos, cs) && return memoized_halos[cs]
    halo_cs = CoordinateSystem{S}(uniquename("halo_" * name(cs)))
    memoized_halos[cs] = halo_cs

    els = cs.elements
    el_meta = element_metadata(cs)
    halo_meta = filter(layer_inclusion(only_layers, ignore_layers), unique(el_meta))

    for meta in halo_meta
        meta_els = els[el_meta .== meta]

        # Get halo of polygons together
        poly::Vector{AbstractPolygon{S}} = filter(x -> x isa AbstractPolygon, meta_els)
        place!.(halo_cs, halo(poly, outer_delta, inner_delta), meta)

        # Get the remaining halos
        isnotpoly = filter(x -> !(x isa AbstractPolygon), meta_els)
        !isempty(isnotpoly) &&
            place!.(halo_cs, halo(isnotpoly, outer_delta, inner_delta), meta)
    end

    # Add halos of references recursively
    refs = map(cs.refs) do ref
        newref = sref(
            halo(
                ref.structure,
                outer_delta,
                inner_delta;
                only_layers=only_layers,
                ignore_layers=ignore_layers,
                memoized_halos=memoized_halos
            ),
            transformation(ref)
        )
        return newref
    end
    halo_cs.refs = refs

    return halo_cs
end

function halo(
    cs::Cell{S},
    outer_delta,
    inner_delta=nothing;
    only_layers=[],
    ignore_layers=[],
    memoized_halos=Dict{GeometryStructure, GeometryStructure}()
) where {S}
    haskey(memoized_halos, cs) && return memoized_halos[cs]
    halo_cs = Cell{S}(uniquename("halo_" * name(cs)))
    memoized_halos[cs] = halo_cs

    els = cs.elements
    el_meta = element_metadata(cs)
    halo_meta = filter(layer_inclusion(only_layers, ignore_layers), unique(el_meta))
    for me in halo_meta
        idx = (cs.element_metadata .== me)
        meta_els = els[idx]
        render!.(halo_cs, halo(meta_els, outer_delta, inner_delta), me)
    end

    refs = map(cs.refs) do ref
        newref = sref(
            halo(
                ref.structure,
                outer_delta,
                inner_delta,
                only_layers=only_layers,
                ignore_layers=ignore_layers,
                memoized_halos=memoized_halos
            ),
            transformation(ref)
        )
        return newref
    end
    addref!.(halo_cs, refs)

    return halo_cs
end

end #module
