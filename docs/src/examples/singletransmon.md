# Single transmon with readout resonator

In this example, we demonstrate the layout and solid model generation of a single transmon
coupled to a readout resonator. In particular, this example demonstrates the additional steps
necessary to generate a mesh suitable for electromagnetic simulation using external
applications, such as the [*Palace*](https://github.com/awslabs/palace) solver. This example
constructs the geometry used in [the *Palace* release documentation](https://aws.amazon.com/blogs/quantum-computing/aws-releases-open-source-software-palace-for-cloud-based-electromagnetics-simulations-of-quantum-computing-hardware/) and is aimed at
demonstrating the power of the [`SolidModels`](../schematicdriven/solidmodels.md)
capability.

The full code for this example can be found [in `examples/SingleTransmon/SingleTransmon.jl` in the `DeviceLayout.jl` repository](https://github.com/aws-cqc/DeviceLayout.jl/blob/main/examples/SingleTransmon/SingleTransmon.jl). Components and process technology are drawn from the
[ExamplePDK](examplepdk.md).

## Overview

The single transmon schematic design is simple, consisting of three components: an
[`ExampleRectangleTransmon`](./examplepdk.md#Transmons), an
[`ExampleClawedMeanderReadout`](./examplepdk.md#ReadoutResonators), and a coplanar waveguide path. The path additionally demonstrates how to terminate a path using lumped ports.

Once the mesh has been generated, `configfile(sm::SolidModel; palace_build=nothing)` will
ingest the `SolidModel` and return a dictionary representing a *Palace* configuration. If
the optional `palace_build` directory is also passed then the configuration file will be
validated using the *Palace* provided JSON schema.

The generated configuration and mesh files can then be run directly from within `julia`
using the `palace_job(config::Dict; palace_build, np=0, nt=1)` which will write the `config`
to `config.json` and if `np > 0` attempt to run the *Palace* job from within the `julia`
shell mode. If `np>0` the configuration file will be written to disk ready for a manual call
to *Palace* outside of `julia`.

These three functions are all wrapped together in `main(palace_build, np=0)` which will
build the `SolidModel`, and write the mesh and configuration file to disk, before optionally
attempting to run *Palace*.

## SolidModel construction

To run the example we `include` the `example/SingleTransmon/SingleTransmon.jl` file and
then call `single_transmon()`. This will construct the schematic and render the design to
`SolidModel`. To visualize the resulting model in `Gmsh`, we then call
`SolidModels.gmsh.fltk.run()`.

```julia
using DeviceLayout
include("examples/SingleTransmon/SingleTransmon.jl")
sm = SingleTransmon.single_transmon()
SolidModels.gmsh.fltk.run() # Opens Gmsh GUI
```

Below, we show the mesh for the metal surfaces. You can also (barely) see the lumped ports at the ends of the readout line in red.

```@raw html
<img src="../../assets/single_transmon_mesh.png"/>
```

The function exposes a number of parameters for the definition of the transmon and readout
resonator, demonstrating how the kernel for an automated parameter search might be
established. Aside from the design parameters there are two processing arguments:
`save_mesh::Bool`, which performs the meshing of the `SolidModel` and writes the resulting
mesh to disk as `single_transmon.msh2` in the Gmsh 2.2 mesh format, and
`save_gds::Bool`, which generates a GDS of the schematic and writes the resulting design to
`single_transmon.gds`.

The schematic design is significantly simplied compared to that of the [quantum processor
example](qpu17.md), but has the same general structure: the three components are
added to the schematic graph and attached to each other before being placed by `plan` into a floorplan, which is then furnished with air bridges.

What differs from the processor example is the construction of the `SolidModel`. When we `render!` the floorplan to a `SolidModel`, instead of a `LayoutTarget` we provide a [`SolidModelTarget`](@ref SchematicDrivenLayout.SolidModelTarget) containing information about how to generate the 3D model. The `SolidModelTarget` is based on this definition from `ExamplePDK`:

```julia
"""
    const SINGLECHIP_SOLIDMODEL_TARGET::SolidModelTarget

A `Target` for rendering to a `SolidModel` using the `ExamplePDK`'s process technology.

Contains rendering options and postrendering operations to create a solid model suitable
for simulation of a single-chip device (as opposed to a flipchip device).
"""
const SINGLECHIP_SOLIDMODEL_TARGET = SolidModelTarget(
    EXAMPLE_SINGLECHIP_TECHNOLOGY; # Thickness and height define z-height and extrusions
    simulation=true, # Optional simulation-only geometry entities will be rendered
    bounding_layers=[:simulated_area], # SIMULATED_AREA defines the simulation bounds
    substrate_layers=[:chip_area], # CHIP_AREA will be extruded downward
    indexed_layers=[:port, :lumped_element, :integration], # Automatically index these layers
    postrender_ops=[ # Manual definition of operations to run after 2D rendering
        (   # Get metal ground plane by subtracting negative from writeable area
            "metal", # Output group name
            SolidModels.difference_geom!, # Operation
            ("writeable_area", "metal_negative", 2, 2), # (object, tool, object_dim, tool_dim)
            :remove_object => true # Remove "writeable_area" group after operation
        ),
        (   # Then add any positive back in
            "metal",
            SolidModels.union_geom!,
            ("metal", "metal_positive", 2, 2),
            :remove_tool => true
        ),
        (   # Define a bulk physical group for all the substrates in the domain.
            "substrate",
            SolidModels.union_geom!,
            ("chip_area_extrusion", "chip_area_extrusion", 3, 3),
            :remove_object => true,
            :remove_tool => true
        ),
        (   # Define the vacuum domain as the remainder of the simulation domain.
            "vacuum",
            SolidModels.difference_geom!,
            ("simulated_area_extrusion", "substrate", 3, 3)
        ),
        # Generate staple bridges in "bridge_metal" group
        SolidModels.staple_bridge_postrendering(;
            base="bridge_base",
            bridge="bridge",
            bridge_height=10μm # Exaggerated, for visualization
        )...,
        (   # Union of all physical metal
            "metal",
            SolidModels.union_geom!,
            ("metal", "bridge_metal"),
            :remove_object => true,
            :remove_tool => true
        ),
        (("metal", SolidModels.difference_geom!, ("metal", "port")))
    ],
    # We only want to retain physical groups that we will need for specifying boundary
    # conditions in the physical domain.
    retained_physical_groups=[
        ("vacuum", 3),
        ("substrate", 3),
        ("metal", 2),
        ("exterior_boundary", 2)
    ]
)
```

We modify it slightly to retain `port_1`, `port_2`, and `lumped_element` physical groups in order to use those in the simulation configuration.

## Configuration

Next, we generate a dictionary defining a *Palace* configuration, specifying the problem type, model, materials, boundary conditions, and solver settings. We can also validate the configuration if we have a path to a *Palace* build, which contains the schema for validation. For more details on configuration, see [the *Palace* documentation](https://awslabs.github.io/palace/stable/config/config/).

For this example, we define most of the configuration by hand. However, we need to identify which materials and boundary conditions apply to which volumes and surfaces in the model.

By construction, the `SolidModel` contains physical groups identifying the volumes and surfaces we're interested in. Using `SolidModels.attributes`, we can get a dictionary mapping the names of the physical groups to the integer "attribute" identifying the corresponding entities in the mesh. We can then use this dictionary to populate the configuration automatically according to the physical intent behind the model, without having to worry about whether the underlying geometry or attribute numbering might change.

```julia
"""
    configfile(sm::SolidModel; palace_build=nothing)

Given a `SolidModel`, assemble a dictionary defining a configuration file for use within
Palace.

  - `sm`: The `SolidModel`from which to construct the configuration file
  - `palace_build = nothing`: Path to a Palace build directory, used to perform validation of
    the configuration file. If not present, no validation is performed.
"""
function configfile(sm::SolidModel; palace_build=nothing)
    attributes = SolidModels.attributes(sm)

    config = Dict(
        "Problem" => Dict(
            "Type" => "Eigenmode",
            "Verbose" => 2,
            "Output" => joinpath(@__DIR__, "postpro/single-transmon")
        ),
        "Model" => Dict(
            "Mesh" => joinpath(@__DIR__, "single_transmon.msh2"),
            "L0" => DeviceLayout.ustrip(m, 1SolidModels.STP_UNIT), # um is Palace default; record it anyway
            "Refinement" => Dict(
                "MaxIts" => 0 # Increase to enable AMR
            )
        ),
        "Domains" => Dict(
            "Materials" => [
                Dict(
                    # Vaccuum
                    "Attributes" => [attributes["vacuum"]],
                    "Permeability" => 1.0,
                    "Permittivity" => 1.0
                ),
                Dict(
                    # Sapphire
                    "Attributes" => [attributes["substrate"]],
                    "Permeability" => [0.99999975, 0.99999975, 0.99999979],
                    "Permittivity" => [9.3, 9.3, 11.5],
                    "LossTan" => [3.0e-5, 3.0e-5, 8.6e-5],
                    "MaterialAxes" =>
                        [[0.8, 0.6, 0.0], [-0.6, 0.8, 0.0], [0.0, 0.0, 1.0]]
                )
            ],
            "Postprocessing" => Dict(
                "Energy" => [Dict("Index" => 1, "Attributes" => [attributes["substrate"]])]
            )
        ),
        "Boundaries" => Dict(
            "PEC" => Dict("Attributes" => [attributes["metal"]]),
            "Absorbing" => Dict(
                "Attributes" => [attributes["exterior_boundary"]],
                "Order" => 1
            ),
            "LumpedPort" => [
                Dict(
                    "Index" => 1,
                    "Attributes" => [attributes["port_1"]],
                    "R" => 50,
                    "Direction" => "+X"
                ),
                Dict(
                    "Index" => 2,
                    "Attributes" => [attributes["port_2"]],
                    "R" => 50,
                    "Direction" => "+X"
                ),
                Dict(
                    "Index" => 3,
                    "Attributes" => [attributes["lumped_element"]],
                    "L" => 14.860e-9,
                    "C" => 5.5e-15,
                    "Direction" => "+Y"
                )
            ]
        ),
        "Solver" => Dict(
            "Order" => 1,
            "Eigenmode" => Dict("N" => 2, "Tol" => 1.0e-6, "Target" => 2, "Save" => 2),
            "Linear" => Dict("Type" => "Default", "Tol" => 1.0e-7, "MaxIts" => 500)
        )
    )

    if !isnothing(palace_build)
        # Load the json schema and validate the configuration
        schema_dir = joinpath(palace_build, "bin", "schema")
        schema = Schema(
            JSON.parsefile(joinpath(schema_dir, "config-schema.json"));
            parent_dir=schema_dir
        )
        validate(schema, config)
    end

    return config
end
```

The output should look like this:

```jl
Configuration: 0.000082 seconds (327 allocations: 31.391 KiB)
Dict{String, Dict{String, Any}} with 5 entries:
  "Problem"    => Dict("Verbose"=>...)
  "Boundaries" => Dict("LumpedPort"=>...)
  "Model"      => Dict("Refinement"=>...)
  "Domains"    => Dict("Postprocessing"=>...)
  "Solver"     => Dict("Eigenmode"=>...)
```

## *Palace*

Finally, we write the configuration to a file, call *Palace*, and parse the computed eigenfrequencies from its output:

```julia
"""
    palace_job(config::Dict; palace_build, np=0, nt=1)

Given a configuration dictionary, write and optionally run the Palace simulation.

Writes `config.json` to `@__DIR__` in order to pass the configuration into Palace.

  - `config` - A configuration file defining the required fields for a Palace configuration file
  - `palace_build` - Path to a Palace build.
  - `np = 0` - Number of MPI processes to use in the call to Palace. If greater than 0 attempts
    to call palace from within the Julia shell. Requires correct specification of `ENV[PATH]`.
  - `nt = 1` - Number of OpenMp threads to use in the call to Palace (requires Palace built with
    OpenMp)
"""
function palace_job(config::Dict; palace_build, np=0, nt=1)
    # Write the configuration file to json, ready for Palace ingestion
    println("Writing configuration file to $(joinpath(@__DIR__, "config.json"))")
    open(joinpath(@__DIR__, "config.json"), "w") do f
        return JSON.print(f, config)
    end

    if np > 0
        # Call Palace using the generated configuration file.
        # Record the terminal output and any error to files.
        println("Running Palace: stdout sent to log.out, stderr sent to err.out")
        withenv("PATH" => "$(ENV["PATH"]):$palace_build/bin") do
            return run(
                pipeline(
                    ignorestatus(
                        `palace -np $np -nt $nt $(joinpath(@__DIR__,"config.json"))`
                    ),
                    stdout=joinpath(@__DIR__, "log.out"),
                    stderr=joinpath(@__DIR__, "err.out")
                )
            )
        end
        println("Complete.")

        # Extract the computed frequencies
        postprodir = joinpath(@__DIR__, config["Problem"]["Output"])
        freq = CSV.File(joinpath(postprodir, "eig.csv"); header=1) |> DataFrame

        println("Eigenmode Frequencies (GHz): ", freq[:, 2])
    end
    return nothing
end
```

The output should look like this:

```jl
Writing configuration file to /path/to/your/config.json
Running Palace: stdout sent to log.out, stderr sent to err.out
Complete.
Eigenmode Frequencies (GHz): [3.161773657, 4.878135762]
Palace: 52.914418 seconds (669 allocations: 46.422 KiB)
```
