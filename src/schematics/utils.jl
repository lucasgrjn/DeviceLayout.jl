function annotated_geometry(comp::AbstractComponent{T}) where {T}
    cs = geometry(comp)
    for (s, h) in pairs(hooks(comp))
        if h isa Vector
            for (i, hh) in enumerate(h)
                arrow = ArrowAnnotation{T}(; text="$(s)_$i")
                arrow_cs = geometry(arrow)
                f = transformation(hh, hooks(arrow).nock)
                push!(cs.refs, sref(arrow_cs, f))
            end
            continue
        end
        arrow = ArrowAnnotation{T}(; text="$s")
        arrow_cs = geometry(arrow)
        f = transformation(h, hooks(arrow).nock)
        push!(cs.refs, sref(arrow_cs, f))
    end
    return cs
end
