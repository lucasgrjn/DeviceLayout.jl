"""
    ExampleTappedHairpin <: Component

A "hairpin" component with a tap and coupling capacitor before the bend and a meander after.

This component is intended for use in demonstrations.

Diagram shows positive metal path for simplicity, but the path is drawn in the `METAL_NEGATIVE` layer.
Hooks are marked with ⋆ and an arrow in their inward direction.

                              :tap
                               ↓
                               ⋆
                               ↕ tap_cap_coupling_distance
                             █████ ↑                         ⤒
                             █████ tap_cap_length            │
                             █████ ↓                      tap_depth
                              ███  ↕ tap_cap_taper_length    │
         ←————tap_position————→█                             ↓
    :p0→⋆█████████████████████████████████                   —
         ←————————straight_length————————→██
        |← ...                             ██
                                          ██
               ███████████████████████████
              |← ... total path length + assumed_extra_length = total_effective_length

# Parameters

Note that when this component is created within an `ExampleFilteredHairpinReadout`,
most of these defaults are overridden.

  - `name = "res"`: Name of component

  - `total_effective_length = 4mm`: Total effective length, such that the resonant frequency
    occurs when this length is a quarter wavelength at that frequency (assuming some effective index)
  - `assumed_extra_length = 1.0mm`: Assumed extra effective length such that the physical path
    length plus `assumed_extra_length` is `total_effective_length`. For example, bends, bridges,
    and coupling capacitors may affect the effective length.
  - `r_bend = 50μm`: Radius of hairpin bend
  - `style = Paths.CPW(10μm, 10μm)`: Path style
  - `initial_snake = Point(0μm, 0μm)`: If nonzero, add an initial s-curve to this point (in the
    hairpin coordinate system, where the path starts from the origin pointing along the positive x-axis)
  - `straight_length = 1.25μm`: Length of first long straight section in hairpin
  - `bridge = nothing`: `CoordinateSystem` holding the air bridge geometry
  - Parameters for tap and coupling capacitor to other resonator

      + `tap_position = 0.7mm`: Position of tap along initial straight segment
      + `tap_location = -1`: Side of hairpin for tap (+1 for right-hand side starting from qubit)
      + `tap_depth = 35μm`: Distance from center of hairpin CPW to end of coupling capacitor
      + `tap_style = Paths.CPW(5μm, 25μm)`: Style of tap path
      + `tap_cap_taper_length = 5μm`: Length of taper from `tap_style` to `tap_cap_style`
      + `tap_cap_length = 10μm`: Length of capacitive pad after taper
      + `tap_cap_style = Paths.CPW(25μm, 15μm)`: Width and gap-width of coupling capacitor as a CPW style
      + `tap_cap_termination_gap = 5μm`: Gap between coupling capacitor metal and ground in the direction of coupling
      + `tap_cap_coupling_distance = 7.5μm`: Distance from end of capacitor metal to `:tap` hook

# Hooks

  - `p0`: Start of hairpin path
  - `tap`: Distance `tap_cap_coupling_distance` away from the end of the tap capacitor metal,
    with `in_direction` pointing back towards the hairpin
"""
@compdef struct ExampleTappedHairpin <: Component
    name = "res"
    total_effective_length = 4mm
    assumed_extra_length = 1.0mm
    r_bend = 50μm
    style = Paths.CPW(10μm, 10μm)
    initial_snake = Point(0μm, 0μm)
    straight_length = 1.25μm
    bridge = nothing
    # Parameters for tap and coupling capacitor to other resonator
    tap_position = 0.7mm
    tap_location = -1
    tap_depth = 35μm
    tap_style = Paths.CPW(5μm, 25μm)
    tap_cap_taper_length = 5μm
    tap_cap_length = 10μm
    tap_cap_style = Paths.CPW(25μm, 15μm)
    tap_cap_termination_gap = 5μm
    tap_cap_coupling_distance = 7.5μm
end

function SchematicDrivenLayout._geometry!(cs::CoordinateSystem, mr::ExampleTappedHairpin)
    paths = _paths(mr)
    place!.(cs, paths)
    return cs
end

function SchematicDrivenLayout.hooks(mr::ExampleTappedHairpin)
    pa, tap = _paths(mr)
    return (; p0=p0_hook(pa), tap=p1_hook(tap))
end

# Helper function for generating paths
function _paths(mr::ExampleTappedHairpin)
    pa = Path(nm; metadata=METAL_NEGATIVE, name=uniquename("hairpin"))
    !iszero(mr.initial_snake) && # If applicable, make an s-curve to `initial_snake`
        route!(
            pa,
            mr.initial_snake,
            0°,
            Paths.StraightAnd45(min_bend_radius=mr.r_bend),
            mr.style
        )

    # First straight leg of hairpin
    straight!(pa, mr.tap_position - Paths.extent(mr.tap_style), mr.style)
    tap = tap!(pa, mr.tap_style; location=mr.tap_location)
    straight!(pa, mr.straight_length - mr.tap_position - Paths.extent(mr.tap_style))
    # Bend and meander with up to straight_length until the total effective length is met
    turn!(pa, sign(mr.tap_location) * 180°, mr.r_bend)
    remainder = mr.total_effective_length - mr.assumed_extra_length - pathlength(pa)
    meander!(pa, remainder, mr.straight_length, mr.r_bend, -sign(mr.tap_location) * 180°)
    terminate!(pa, gap=0μm, rounding=(mr.style.gap / 2)) # Round the ends of the short
    add_bridges!(pa, mr.bridge)
    # Extend `tap` to make a coupling capacitor
    straight!(
        tap,
        mr.tap_depth - Paths.extent(mr.style) - mr.tap_cap_length - mr.tap_cap_taper_length
    )
    straight!(tap, mr.tap_cap_taper_length, Paths.Taper())
    straight!(tap, mr.tap_cap_length, mr.tap_cap_style)
    terminate!(tap; gap=mr.tap_cap_termination_gap)
    # Add NoRender segment to get to the hook at coupling_distance away from the metal termination
    straight!(
        tap,
        mr.tap_cap_coupling_distance - mr.tap_cap_termination_gap,
        Paths.NoRender()
    )

    return pa, tap
end
