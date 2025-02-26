using Test, BenchmarkTools
const SUITE = BenchmarkGroup()

using DeviceLayout, Pkg, FileIO, Unitful, DeviceLayout.PreferredUnits
const um = μm

include(joinpath(@__DIR__, "clipping.jl"))
