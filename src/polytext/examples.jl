const scripted_equation = "H=ħωσ_x+Jσ_y

c^2=a^2+b^2

|Ψ>=a^†a|ψ_1>+b^†b|ψ_2>

e^{𝚤π}=-1 and c_{vac}≈3x10^8 m/s"

const test_string = string(
    "Keyboard characters:

    ABCDEFGHIJKLMNOPQRSTUVWXYZ
    abcdefghijklmnopqrstuvwxyz

    ~!@#\$%^&*()_+
    `1234567890-=

    []\\{}|;':\",./<>?

    Non-keyboard characters:

    αβγδϵηθκλμνπρστϕχψωħ
    ΩΣΞΓΠΨΔΛΘΦ
    𝚤∞÷√±∓≠≡≈⟂∠°∂∇†□░█\n\n",
    length(lcd),
    " characters supported so far.

Line limit demonstration (newline character not present):
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Aliquam in enim vestibulum, laoreet ligula at, convallis nisl. Curabitur elit mi, luctus a semper sed, euismod sed turpis. Nunc ac arcu egestas, tristique leo vitae, pellentesque augue. Vivamus massa urna, varius quis scelerisque ac, imperdiet non magna. Curabitur id rhoncus nisl. Cras consequat vulputate mauris, sit amet congue odio. Sed posuere ullamcorper libero, id efficitur diam auctor quis. Morbi ac neque lectus. Maecenas ultrices placerat justo, id sollicitudin velit dapibus laoreet. Mauris sodales consectetur mi eget suscipit. Morbi eu rutrum turpis. In sed dolor eu purus venenatis feugiat. Maecenas lacinia dui vel consequat venenatis. Aenean viverra, quam nec tempus iaculis, velit libero laoreet ligula, id hendrerit velit lorem at velit."
)

const reference_test_string = "aaaa
bbbb
cccc
dddd
eeee
ffff
gggg
♇♇♇♇" # End with invalid character to test catch fallback

"""
    scripted_demo(save_path = joinpath(homedir(),"Desktop","scripted.gds"), flatten = false)

Demo script for demonstrating the use of the `scripting` parameter in `polytext!`.
`flatten` can flatten the cells before saving (for SVG output).
"""
function scripted_demo(
    save_path=joinpath(homedir(), "Desktop", "scripted.gds"),
    flatten=false
)
    reset_uniquename!()
    c = Cell("scripted", nm)
    sty = DotMatrix(; pixelsize=1μm)
    polytext!(c, scripted_equation, sty; scripting=true)
    flatten && flatten!(c)
    return save(save_path, c)
end

"""
    characters_demo(save_path = joinpath(homedir(),"Desktop","characters.gds"), flatten = false)

Demo script for demonstrating the available characters in `polytext!` and the `linelimit`
parameter in use. `flatten` can flatten the cells before saving (for SVG output).
"""
function characters_demo(
    save_path=joinpath(homedir(), "Desktop", "characters.gds"),
    flatten=false
)
    reset_uniquename!()
    c = Cell("characters", nm)
    sty = DotMatrix(; pixelsize=1μm)
    polytext!(c, test_string, sty; linelimit=80)
    flatten && flatten!(c)
    return save(save_path, c)
end

"""
    referenced_characters_demo(save_path = joinpath(homedir(),"Desktop","referenced_characters.gds");
        verbose_override = false)

Demo script for demonstrating the memory saving ability of keeping CellReferences for
previously used characters in `polytext!`. Nothing is printed if `verbose_override` is `true`.
"""
function referenced_characters_demo(
    save_path=joinpath(homedir(), "Desktop", "referenced_characters.gds");
    verbose_override=false
)
    reset_uniquename!()
    c = Cell("referenced_characters", nm)
    sty = DotMatrix(; pixelsize=1μm)
    polytext!(c, reference_test_string, sty; verbose=(!verbose_override))
    return save(save_path, c)
end
