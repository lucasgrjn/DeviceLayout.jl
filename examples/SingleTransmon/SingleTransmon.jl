
module SingleTransmon
# Using a module to let us cleanly include file in test script without worrying about namespaces

using FileIO, CSV, DataFrames, JSON, JSONSchema
using DeviceLayout, DeviceLayout.SchematicDrivenLayout, DeviceLayout.PreferredUnits
import .SchematicDrivenLayout.ExamplePDK
import .SchematicDrivenLayout.ExamplePDK: LayerVocabulary, L1_TARGET, add_bridges!
using .ExamplePDK.Transmons, .ExamplePDK.ReadoutResonators
import .ExamplePDK.SimpleJunctions: ExampleSimpleJunction

"""
    single_transmon(
        w_shield=2μm,
        claw_gap=6μm,
        w_claw=34μm,
        l_claw=121μm,
        cap_width=24μm,
        cap_length=620μm,
        cap_gap=30μm,
        n_meander_turns=5,
        hanger_length=500μm,
        bend_radius=50μm,
        save_mesh::Bool=false,
        save_gds::Bool=false)

Generate a SolidModel and mesh for a single transmon design, using a rectangular transmon island and claw resonator.
"""
function single_transmon(;
    w_shield=2μm,
    claw_gap=6μm,
    w_claw=34μm,
    l_claw=121μm,
    cap_width=24μm,
    cap_length=620μm,
    cap_gap=30μm,
    n_meander_turns=5,
    hanger_length=500μm,
    bend_radius=50μm,
    save_mesh::Bool=false,
    save_gds::Bool=false
)
    #### Reset name counter for consistency within a Julia session
    reset_uniquename!()
    # Compute implicit parameters
    w_grasp = cap_width + 2 * cap_gap
    total_height = 444μm + hanger_length + (3 + n_meander_turns * 2) * bend_radius

    #### Create abstract components
    cpw_width = 10μm
    cpw_gap = 6μm
    PATH_STYLE = Paths.SimpleCPW(cpw_width, cpw_gap)
    BRIDGE_STYLE = ExamplePDK.bridge_geometry(PATH_STYLE)

    qubit = ExampleRectangleTransmon(;
        jj_template=ExampleSimpleJunction(),
        name="qubit",
        cap_length,
        cap_gap,
        cap_width
    )
    rres = ExampleClawedMeanderReadout(;
        name="rres",
        coupling_length=400μm,
        w_shield,
        w_claw,
        l_claw,
        claw_gap,
        w_grasp,
        n_meander_turns,
        total_height,
        hanger_length,
        bend_radius,
        bridge=BRIDGE_STYLE
    )

    ##### Build schematic graph
    g = SchematicGraph("single-transmon")
    qubit_node = add_node!(g, qubit)
    rres_node = fuse!(g, qubit_node, rres)
    # Equivalent to `fuse!(g, qubit_node=>:readout, rres=>:qubit)`
    # because `matching_hooks` was implemented for that component pair

    ## Readout path
    readout_length = 2700μm
    p_readout = Path(
        Point(0μm, 0μm);
        α0=π / 2,
        name="p_ro",
        metadata=LayerVocabulary.METAL_NEGATIVE
    )
    straight!(p_readout, readout_length / 2, PATH_STYLE)
    straight!(p_readout, readout_length / 2, PATH_STYLE)

    ## Readout lumped ports - one at each end
    csport = CoordinateSystem(uniquename("port"), nm)
    render!(
        csport,
        only_simulated(centered(Rectangle(cpw_width, cpw_width))),
        LayerVocabulary.PORT
    )
    attach!(p_readout, sref(csport), cpw_width, i=1) # @ start
    attach!(p_readout, sref(csport), readout_length / 2 - cpw_width, i=2) # @ end
    p_readout_node = add_node!(g, p_readout)
    ## Attach resonator to feedline
    # Instead of `fuse!` we use a schematic-based `attach!` method to place it along the path
    # Syntax is a mix of `fuse!` and how we attached the ports above
    attach!(g, p_readout_node, rres_node => :feedline, 0mm, location=1)

    ## Position the components.
    floorplan = plan(g)
    add_bridges!(floorplan, BRIDGE_STYLE, spacing=300μm) # Add bridges to paths

    #### Prepare solid model
    # Specify the extent of the simulation domain.
    substrate_x = 4mm
    substrate_y = 3.7mm

    center_xyz = DeviceLayout.center(floorplan)
    chip = centered(Rectangle(substrate_x, substrate_y), on_pt=center_xyz)
    sim_area = centered(Rectangle(substrate_x, substrate_y), on_pt=center_xyz)

    render!(floorplan.coordinate_system, sim_area, LayerVocabulary.SIMULATED_AREA)
    render!(floorplan.coordinate_system, sim_area, LayerVocabulary.WRITEABLE_AREA)
    render!(floorplan.coordinate_system, chip, LayerVocabulary.CHIP_AREA)

    check!(floorplan)

    # Need to pass generated physical group names so they can be retained
    tech = ExamplePDK.singlechip_solidmodel_target("port_1", "port_2", "lumped_element")
    sm = SolidModel("test", overwrite=true)

    # Adjust mesh_scale to increase the resolution of the mesh, < 1 will result in greater
    # resolution near edges of the geometry.
    meshing_parameters = SolidModels.MeshingParameters(
        mesh_scale=1.0,
        mesh_order=2,
        options=Dict("General.Verbosity" => 1) # General Gmsh option input
    )
    render!(sm, floorplan, tech, meshing_parameters=meshing_parameters)

    if save_mesh
        SolidModels.gmsh.model.mesh.generate(3) # runs without error
        save(joinpath(@__DIR__, "single_transmon.msh2"), sm)
    end

    if save_gds
        # Render to GDS as well.
        c = Cell("single_transmon", nm)
        render!(c, floorplan, L1_TARGET, strict=:no)
        flatten!(c)
        save(joinpath(@__DIR__, "single_transmon.gds"), c)
    end
    return sm
end

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

"""
    main(palace_build, np=0)

Given a build of Palace found at `palace_build`, assemble a `SolidModel` of the single
transmon example, mesh the geometry, define a configuration file for the geometry, and if
`np > 0`, launch Palace from the provided `palace_build`.`
"""
function main(palace_build, np=0)
    # Construct the SolidModel
    @time "SolidModel + Meshing" sm = single_transmon(save_mesh=true)
    # Assemble the configuration
    @time "Configuration" config = configfile(sm; palace_build)
    # Call Palace
    @time "Palace" palace_job(config; palace_build, np)
end

end # module
