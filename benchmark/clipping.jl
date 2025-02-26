SUITE["clipping"] = BenchmarkGroup()
SUITE["clipping"]["difference2d"] = BenchmarkGroup()

for name in [
    :difference2d_100_square,
    :difference2d_961_square,
    :difference2d_3600_square,
    :difference2d_10000_square,
    :difference2d_19881_square
]
    @eval str = string($(QuoteNode(name)))
    @eval $name = load(joinpath(dirname(@__FILE__), "data", str * ".gds"))["main"]
    @eval b = elements($name)
    @eval c = element_metadata($name)
    @eval a = splice!(b, findall(m -> m == GDSMeta(0, 0), c))

    names = split(str, "_")
    n1 = names[1]
    n2 = join(names[2:end], "_")
    SUITE["clipping"][n1][n2] = @benchmarkable difference2d($(a), $(b))
end

difference2d_8000_skew =
    load(joinpath(dirname(@__FILE__), "data", "difference2d_8000_skew.gds"))["main"]
b = elements(difference2d_8000_skew)
c = element_metadata(difference2d_8000_skew)
a = splice!(b, findall(m -> m == GDSMeta(0, 0), c))

SUITE["clipping"]["difference2d"]["8000_skew"] = @benchmarkable difference2d($(a), $(b))
