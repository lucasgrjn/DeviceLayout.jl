using PrecompileTools

@static if unit_preference == "NoUnits"
    @setup_workload begin
        @compile_workload begin
            cs = CoordinateSystem("test")
            r = Polygons.Rounded(simple_tee(1, 10), 1)
            place!(cs, r, :test)
        end
    end
else
    @setup_workload begin
        @compile_workload begin
            cs = CoordinateSystem("test", nm)
            cs2 = CoordinateSystem("attachment", nm)
            r = Polygons.Rounded(simple_tee(1μm, 10μm), 1μm)
            place!(cs2, r, :test)

            pa = Path(nm)
            straight!(pa, 100μm, Paths.SimpleCPW(10μm, 6μm))
            turn!(pa, 45°, 100μm)
            attach!(pa, sref(cs2), pathlength(pa[end]))
            pa.metadata = SemanticMeta(:test)

            addref!(cs, sref(pa, rot=45°))
            c = Cell(cs, map_meta=_ -> GDSMeta())
        end
    end
end
