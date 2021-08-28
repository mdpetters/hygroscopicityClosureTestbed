using Cairo
using Fontconfig
using Gadfly
using DataFrames
using Colors
using Printf

include("AIOMFAC.jl")
using .AIOMFAC

Tk = 298.0
components = (("Na", 1), ("Cl", 1))
writeinput("input_0001.txt", Tk, components)
run(`./AIOMFAC Inputfiles/input_0001.txt`)
aw, k = parseoutput("AIOMFAC_output_0001.txt")
rhoNaCl = 2.16
df1 = DataFrame(aw=aw, k=k .* rhoNaCl, comp="NaCl")

components = (("H", 2), ("SO4", 1))
writeinput("input_0001.txt", Tk, components)
run(`./AIOMFAC Inputfiles/input_0001.txt`)
aw, k = parseoutput("AIOMFAC_output_0001.txt")
rhoH2SO4 = 1.83
df2 = DataFrame(aw=aw, k=k .* rhoH2SO4, comp="H2SO4")

components = (("Na", 2), ("SO4", 1))
writeinput("input_0001.txt", Tk, components)
run(`./AIOMFAC Inputfiles/input_0001.txt`)
aw, k = parseoutput("AIOMFAC_output_0001.txt")
rhoSodSul = 2.66
df3 = DataFrame(aw=aw, k=k .* rhoSodSul, comp="Na2SO4")

components = (("NH4", 2), ("SO4", 1))
writeinput("input_0001.txt", Tk, components)
run(`./AIOMFAC Inputfiles/input_0001.txt`)
aw, k = parseoutput("AIOMFAC_output_0001.txt")
rhoAS = 1.76
df4 = DataFrame(aw=aw, k=k .* rhoAS, comp="(NH4)2SO4")

components = (("NH4", 1), ("HSO4", 1))
writeinput("input_0001.txt", Tk, components)
run(`./AIOMFAC Inputfiles/input_0001.txt`)
aw, k = parseoutput("AIOMFAC_output_0001.txt")
rhoABS = 1.78
df5 = DataFrame(aw=aw, k=k .* rhoABS, comp="NH4HSO4")

df = [df1;df2;df3;df4;df5]
colors = ["black", "darkred", "steelblue3", "darkgoldenrod", "grey"]
ylabels = [0, 1, 2, 3]
p = plot(df, x=:aw, y=:k, color=:comp, Geom.line,
	Guide.xlabel("Water activity (-)"),
	Guide.ylabel("κᵥ"),
	Guide.colorkey(title="Composition"),
	Guide.xticks(ticks=0.4:0.1:1.0),
	Guide.yticks(ticks=0.0:0.2:3.0),
    Scale.y_continuous(labels = y -> y in ylabels ? @sprintf("%.1f", y) : ""),
	Scale.color_discrete_manual(colors...),
	Coord.cartesian(xmin=0.4, xmax=1, ymin=0.0, ymax=3)
)

set_default_plot_size(4.5inch, 3inch)
draw(PNG("validate.png", dpi=300), p)

