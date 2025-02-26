# DeviceLayout.jl

DeviceLayout.jl is a [Julia](http://julialang.org) package for computer-aided design (CAD) of quantum integrated circuits developed at the AWS Center for Quantum Computing. The package supports the generation of 2D layouts and 3D models of complex devices using a low-level geometry interface together with a high-level schematic-driven workflow.

## Installation

Julia can be downloaded [here](https://julialang.org/downloads/). We support Julia v1.9 or later.

From Julia, install DeviceLayout.jl using the built-in package manager:

```julia
import Pkg
Pkg.activate(".") # Activates an environment in the current directory
Pkg.add("DeviceLayout")
```

We recommend [using an environment for each project](https://julialang.github.io/Pkg.jl/v1/environments/) rather than installing packages in the default environment.

## Quick start

Let's mock up a transmission line with two launchers and some bridges across the
transmission line. We begin by making a cell with a rectangle in it:

```@example 1
using DeviceLayout, DeviceLayout.PreferredUnits, FileIO

cr = Cell("rect", nm)
r = centered(Rectangle(20μm, 40μm))
render!(cr, r, GDSMeta(1, 0))
save("units_rectonly.svg", cr; layercolors=Dict(0 => (0, 0, 0, 1), 1 => (1, 0, 0, 1)));
nothing; # hide
```

```@raw html
<img src="units_rectonly.svg" style="width:1in;"/>
```

A rectangle made with a width and height parameter will default to having its lower-left
corner at the origin. [`centered`](@ref) will return a rectangle that is centered about the origin
instead.

The rectangle is then rendered into the cell. `GDSMeta(1)` indicates the target GDSII layer. You
can also specify the GDSII datatype as a second argument, e.g. `GDSMeta(1,0)`.

In another cell, we make the transmission line with some launchers on both ends:

```@example 1
p = Path(μm)
sty = launch!(p)
straight!(p, 500μm, sty)
turn!(p, π / 2, 150μm)
straight!(p, 500μm)
launch!(p)
cp = Cell("pathonly", nm)
render!(cp, p, GDSMeta(0))
save("units_pathonly.svg", cp; layercolors=Dict(0 => (0, 0, 0, 1), 1 => (1, 0, 0, 1)));
nothing; # hide
```

```@raw html
<img src="units_pathonly.svg" style="width: 3in;"/>
```

Finally, let's put bridges across the feedline:

```@example 1
turnidx = Int((length(p) + 1) / 2) - 1 # the first straight segment of the path
simplify!(p, turnidx .+ (0:2))
attach!(
    p,
    CellReference(cr, Point(0.0μm, 0.0μm)),
    (40μm):(40μm):((pathlength(p[turnidx])) - 40μm),
    i=turnidx
)
c = Cell("decoratedpath", nm)
render!(c, p, GDSMeta(0))
save("units.svg", flatten(c); layercolors=Dict(0 => (0, 0, 0, 1), 1 => (1, 0, 0, 1)));
nothing; # hide
```

```@raw html
<img src="units.svg" style="width: 3in;"/>
```

You can save a `Cell` to a GDS file for lithography or an SVG for vector graphics by using
`save` with an appropriate extension:

```julia
save("/path/to/myoutput.gds", c)
save("/path/to/myoutput.svg", c)
```

SVG support is experimental, but it is used in generating the graphics you see in this documentation. If you use [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/#The-Julia-REPL) (read-eval-print loop) provided by the extension [Julia for Visual Studio Code](https://www.julia-vscode.org/),
rendered cells are previewed in a separate tab. If you use Jupyter/IJulia, rendered
cells are automatically returned as a result.

## Performance and workflow tips

[KLayout](https://www.klayout.de/) is a free (GPL v2+) GDS viewer/editor. It watches
its open files for changes, making it easy to use as a fast previewer alongside DeviceLayout.jl.

The recommended IDE for Julia is [Visual Studio Code](https://code.visualstudio.com/) with the [Julia for Visual Studio Code extension](https://www.julia-vscode.org/).

Since Julia has a just-in-time compiler, the first time code is executed may take much
longer than any other times. This means that a lot of time will be wasted repeating
compilations if you run DeviceLayout.jl in a script like you would in other languages. For
readability, it is best to split up your CAD code into functions that have clearly named
inputs and perform a well-defined task.

It is also best to avoid writing statements in global scope. In other words, put most of
your code in a function. Your CAD script should ideally look like the following:

```julia
using DeviceLayout, DeviceLayout.PreferredUnits, FileIO

function subroutine1()
    # render some thing
end

function subroutine2()
    # render some other thing
end

function main()
    # my cad code goes here: do all of the things
    subroutine1()
    subroutine2()
    return save("/path/to/out.gds", ...)
end

main() # execute main() at end of script.
```

In a typical workflow, you'll have a text editor open alongside a Julia REPL. You'll save the above code in a file (e.g., `mycad.jl`) and then run `include("mycad.jl")` from the Julia REPL to generate your pattern.
You'll iteratively revise `mycad.jl` and save your changes.
Subsequent runs should be several times faster than the first, if you `include` the file again from the same Julia session.
