using DeviceLayout, FileIO, Unitful, DeviceLayout.PreferredUnits
const um = μm

for n in (10, 31, 60, 100, 141)
    c = Cell("rect", nm)
    r = Rectangle(8um, 8um)
    render!(c, r, GDSMeta(1, 0))

    x = 10n - 1
    rarr = aref(c, (0:10:x)um, (0:10:x)um)
    main = Cell("main", nm)
    push!(main.refs, rarr)
    flatten!(main)

    rbound = Rectangle(Point(-1, -1)um, Point(x, x)um)
    render!(main, rbound, GDSMeta(0, 0))
    save(joinpath(dirname(@__FILE__), "difference2d_$(n^2)_square.gds"), main)
    # comment above two lines, uncomment below two lines to do differencing locally
    # plgn = difference2d(rbound, elements(main))
    # render!(main, plgn, GDSMeta(2))
end

# difference2d_8000_skew

c = Cell("rect", nm)
r = Rectangle(8um, 8um)
render!(c, r, GDSMeta(1, 0))

rarr = aref(c; dr=Point(1, 10)um, dc=Point(10, 0)um, nr=100, nc=80)
main = Cell("main", nm)
push!(main.refs, rarr)
flatten!(main)

rbound = Rectangle(Point(-1, -1)um, Point(999, 999)um)
render!(main, rbound, GDSMeta(0, 0))
save(joinpath(dirname(@__FILE__), "difference2d_8000_skew.gds"), main)
# comment above two lines, uncomment below two lines to do differencing locally
# plgn = difference2d(rbound, elements(main))
# render!(main, plgn, GDSMeta(2))
