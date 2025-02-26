import DeviceLayout: flushtop, flushleft, flushright, below

import .ExamplePDK: filter_params, tap!
import .ExamplePDK.Transmons: ExampleRectangleTransmon

"""
    struct ExampleClawedMeanderReadout <: Component

Readout resonator consisting of a meander with short and claw-capacitor terminations.

This component is intended for use in demonstrations with `ExampleRectangleTransmon`.

# Parameters

  - `name`: Name of component
  - `style`: Resonator CPW style
  - `total_length`: Total length of resonator
  - `coupling_length`: Length of coupling section
  - `coupling_gap`: Width of ground plane between coupling section and coupled line
  - `bend_radius`: Meander bend radius
  - `n_meander_turns`: Number of meander turns
  - `total_height`: Total height from top hook to bottom hook
  - `hanger_length`: Length of hanger section between coupling section and meander
  - `w_shield`: Width of claw capacitor ground plane shield
  - `w_claw`: Claw trace width
  - `l_claw`: Claw finger length
  - `claw_gap`: Claw capacitor gap
  - `w_grasp`: Width between inner edges of ground plane shield
  - `bridge`: `CoordinateSystem` containing a bridge
"""
@compdef struct ExampleClawedMeanderReadout <: Component
    name            = "rres"
    style           = DeviceLayout.Paths.SimpleCPW(10.0μm, 6.0μm)
    total_length    = 5000μm
    coupling_length = 200μm
    coupling_gap    = 5μm
    bend_radius     = 50μm
    n_meander_turns = 5
    total_height    = 1656μm # from top hook to bottom hook
    hanger_length   = 500μm
    w_shield        = 2μm
    w_claw          = 32μm
    l_claw          = 160μm
    claw_gap        = 6μm
    w_grasp         = 84μm
    bridge          = CoordinateSystem("rresbridge", nm)
end

function SchematicDrivenLayout._geometry!(
    cs::CoordinateSystem,
    rres::ExampleClawedMeanderReadout
)
    params = parameters(rres)
    # Center vertical axis is midpoint of coupling section
    pres = Path(
        Point(
            -params.coupling_length / 2,
            -params.coupling_gap - params.style.gap - params.style.trace / 2
        ),
        α0=0°
    )
    n_bends = 3 + 2 * params.n_meander_turns # number of 90 degree bends
    arm_length = (
        params.total_height - params.hanger_length - n_bends * params.bend_radius -
        params.coupling_gap - params.style.gap - params.style.trace / 2 -
        params.w_shield - 2 * params.claw_gap - params.w_claw
    )

    # Length of straight sections in meander
    straight_length =
        (
            params.total_length - 3 * params.coupling_length / 2 -
            n_bends * pi * params.bend_radius / 2 - arm_length - params.hanger_length
        ) / params.n_meander_turns

    bb_cs = params.bridge
    ###

    straight!(pres, params.coupling_length, params.style)
    turn!(pres, -90°, params.bend_radius)
    straight!(pres, params.hanger_length)
    attach!(pres, CoordinateSystemReference(bb_cs), params.hanger_length / 2)
    turn!(pres, -90°, params.bend_radius)
    # Center of the straight section of meander lines up with coupling midpoint (and claw)
    straight!(pres, straight_length / 2 + params.coupling_length / 2)
    turn!(pres, 180°, params.bend_radius)

    # Start the meander with a full straight section
    meander_length =
        (params.n_meander_turns - 1) * (straight_length + pi * params.bend_radius) +
        straight_length / 2 - params.bend_radius
    meander!(pres, meander_length, straight_length, params.bend_radius, -180°)
    turn!(pres, -90°, params.bend_radius)
    straight!(pres, arm_length)
    attach!(pres, CoordinateSystemReference(bb_cs), arm_length / 2)

    ### Claw
    w_shield = params.w_shield
    w_claw = params.w_claw
    l_claw = params.l_claw
    claw_gap = params.claw_gap
    arm_trace = params.style.trace
    w_grasp = params.w_grasp # qubit cap trace width + 2 gap width

    p0 = p1(pres.nodes[end].seg)

    claw_hole1 = Rectangle(arm_trace, claw_gap) + p0 + Point(-arm_trace / 2, -claw_gap)

    claw_hole2 =
        Rectangle(w_grasp + 2 * w_shield + 4 * claw_gap + 2 * w_claw, w_claw + 2 * claw_gap)
    claw_hole2 = flushtop(claw_hole2, claw_hole1, centered=true)

    claw_hole3 = Rectangle(w_claw + 2 * claw_gap, w_shield + l_claw + claw_gap)
    claw_hole3 = flushleft(below(claw_hole3, claw_hole2), claw_hole2)

    claw_hole4 = flushright(claw_hole3, claw_hole2)

    claw1 = Rectangle(arm_trace, claw_gap)
    claw1 = flushtop(claw1, claw_hole1, centered=true)

    claw2 = Rectangle(w_grasp + 2 * w_shield + 2 * claw_gap + 2 * w_claw, w_claw)
    claw2 = below(claw2, claw1, centered=true)

    claw3 = Rectangle(w_claw, claw_gap + w_shield + l_claw)
    claw3 = flushleft(below(claw3, claw2), claw2)

    claw4 = flushright(claw3, claw2)

    claw = difference2d(
        [claw_hole1, claw_hole2, claw_hole3, claw_hole4],
        [claw1, claw2, claw3, claw4]
    )

    render!.(cs, [pres, claw], METAL_NEGATIVE)
    return cs
end

"""
    hooks(rres::ExampleClawedMeanderReadout)

`Hook`s for attaching a readout resonator claw to a qubit and coupling section to a feedline.

  - `qubit`: The "palm" of the claw on the outside edge of the "shield". Matches
    `(eq::ExampleRectangleTransmon) => :rres`.
  - `feedline`: A distance `coupling_gap` from the edge of the ground plane, vertically aligned
    with the claw.
"""
function SchematicDrivenLayout.hooks(rres::ExampleClawedMeanderReadout)
    params = parameters(rres)
    qubit_hook = PointHook(Point(zero(params.w_claw), -params.total_height), 90°)
    feedline_hook = PointHook(zero(Point{typeof(params.w_claw)}), -90°)
    return (qubit=qubit_hook, feedline=feedline_hook)
end

SchematicDrivenLayout.matching_hooks(
    ::ExampleRectangleTransmon,
    ::ExampleClawedMeanderReadout
) = (:readout, :qubit)
SchematicDrivenLayout.matching_hooks(
    ::ExampleClawedMeanderReadout,
    ::ExampleRectangleTransmon
) = (:qubit, :readout)
