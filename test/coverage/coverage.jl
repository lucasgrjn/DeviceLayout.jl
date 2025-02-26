using LocalCoverage, Pkg

dir = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.develop(path=dir)
cov = generate_coverage("DeviceLayout", run_test=true)
generate_xml(cov)
Pkg.rm("DeviceLayout")
show(IOContext(stdout, :print_gaps => true), cov)
