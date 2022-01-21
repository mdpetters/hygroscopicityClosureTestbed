using Cairo
using Fontconfig
using Gadfly
using DataFrames
using Colors
using Printf
using CSV
using Lazy
using AIOMFAC

function component2kappa(x)
    aw = x.RH ./ 100.0
    Gm = @. 1.0 / (1.0 - x.w)
    k = (@. (1.0 / aw - 1.0) * (Gm - 1.0))

    return aw, k
end

dfx = CSV.read("../data/Wex2010.txt", DataFrame, header = false)
dfx1 = DataFrame(
    aw = dfx[1:11, 1],
    k = dfx[1:11, 2] ./ 2.01,
    kmin = dfx[12:22, 2] ./ 2.01,
    kmax = (dfx[1:11, 2] + (dfx[1:11, 2] - dfx[12:22, 2])) ./ 2.01,
)

Tk = 298.0
fi = [0.005:0.001:0.009; 0.01:0.01:0.99]
path = "/home/petters/Documents/Research/hygroscopicityClosureTestbed/paper/"
components = [(("Na+", 1), ("Cl-", 1))]
AIOMFAC.writeinput("input_0001.txt", Tk, fi, components, 1.0)
run(`./AIOMFAC Inputfiles/input_0001.txt`)
phase, component = AIOMFAC.load_data("AIOMFAC_output_0001.txt", Tk, fi)
aw, k = component2kappa(component[1])
df1 = DataFrame(aw = aw, k = k, comp = "NaCl/seaspray")

components = [(("H+", 2), ("SO4--", 1))]
AIOMFAC.writeinput("input_0001.txt", Tk, fi, components, 1.0)
run(`./AIOMFAC Inputfiles/input_0001.txt`)
phase, component = AIOMFAC.load_data("AIOMFAC_output_0001.txt", Tk, fi)
aw, k = component2kappa(component[1])
df2 = DataFrame(aw = aw, k = k, comp = "H2SO4")

components = [(("Na+", 2), ("SO4--", 1))]
AIOMFAC.writeinput("input_0001.txt", Tk, fi, components, 1.0)
run(`./AIOMFAC Inputfiles/input_0001.txt`)
phase, component = AIOMFAC.load_data("AIOMFAC_output_0001.txt", Tk, fi)
aw, k = component2kappa(component[1])
df3 = DataFrame(aw = aw, k = k, comp = "Na2SO4")

components = [(("NH4+", 2), ("SO4--", 1))]
AIOMFAC.writeinput("input_0001.txt", Tk, fi, components, 1.0)
run(`./AIOMFAC Inputfiles/input_0001.txt`)
phase, component = AIOMFAC.load_data("AIOMFAC_output_0001.txt", Tk, fi)
aw, k = component2kappa(component[1])
df4 = DataFrame(aw = aw, k = k, comp = "(NH4)2SO4")

components = [(("NH4+", 1), ("HSO4-", 1))]
AIOMFAC.writeinput("input_0001.txt", Tk, fi, components, 1.0)
run(`./AIOMFAC Inputfiles/input_0001.txt`)
phase, component = AIOMFAC.load_data("AIOMFAC_output_0001.txt", Tk, fi)
aw, k = component2kappa(component[1])
df5 = DataFrame(aw = aw, k = k, comp = "NH4HSO4")

df = [df1; df2; df3; df4; df5]
colors = ["black", "darkred", "steelblue3", "darkgoldenrod", "grey", "seagreen"]
ylabels = [0, 0.5, 1]
p = plot(
    layer(
        dfx1,
        x = :aw,
        y = :k,
        ymin = :kmin,
        ymax = :kmax,
        Geom.errorbar,
        Geom.point,
        Theme(default_color = "black"),
    ),
    layer(df, x = :aw, y = :k, color = :comp, Geom.line),
    layer(
        x = [0.5, 1.0],
        y = [0.05, 0.05],
        ymin = [0.005, 0.005],
        ymax = [0.1, 0.1],
        color = ["organic aerosol", "organic aerosol"],
        Geom.line,
        Geom.ribbon,
        Theme(alphas = [0.2], lowlight_color = c -> RGBA{Float32}(c.r, c.g, c.b)),
    ),
    Guide.xlabel("Water activity (-)"),
    Guide.ylabel("κₘ"),
    Guide.colorkey(title = "Composition"),
    Guide.xticks(ticks = 0.5:0.1:1.0),
    Guide.yticks(ticks = 0.0:0.1:1.2),
    Scale.y_continuous(labels = y -> y in ylabels ? @sprintf("%.1f", y) : ""),
    Scale.color_discrete_manual(colors...),
    Theme(plot_padding = [0pt, 11pt, 1pt, 0pt]),
    Coord.cartesian(xmin = 0.5, xmax = 1, ymin = 0.0, ymax = 1.2),
)

set_default_plot_size(4.75inch, 3inch)
draw(PNG("f02.png", dpi = 600), p)
