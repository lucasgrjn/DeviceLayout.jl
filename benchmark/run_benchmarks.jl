using PkgBenchmark
using UnicodePlots, Statistics

include(joinpath(dirname(@__FILE__), "data", "generate_data.jl"))
bench = benchmarkpkg("DeviceLayout")
export_markdown("benchmark.md", bench)
export_markdown(stdout, bench)

diff2dgroup = bench.benchmarkgroup["clipping"]["difference2d"]

for series in ("square",)
    N = Vector{Float64}(undef, 0)
    t = Vector{Float64}(undef, 0)
    for k in keys(diff2dgroup)
        n, v = split(k, '_')
        if v == series
            push!(N, parse(Float64, n))
            push!(t, median(diff2dgroup[k].times) / 1e9)
        end
    end
    idx = sortperm(N)
    display(lineplot(N[idx], t[idx]))
end
