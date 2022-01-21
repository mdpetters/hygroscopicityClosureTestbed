using Cairo
using Fontconfig
using Gadfly
using DataFrames
using MLStyle
using Chain
using Printf
using Colors
using NumericIO

mutable struct SizeDistribution
    A::Any                        # Input parameters [[N1,Dg1,σg1], ...] or DMA
    De::Vector{<:AbstractFloat}   # bin edges
    Dp::Vector{<:AbstractFloat}   # bin midpoints
    ΔlnD::Vector{<:AbstractFloat} # ΔlnD of the grid
    S::Vector{<:AbstractFloat}    # spectral density
    N::Vector{<:AbstractFloat}    # number concentration per bin
    form::Symbol                  # form of the size distribution [:lognormal, ....]
end

function logndist(A, De)
    md(A, x) = @. A[1] / (√(2π) * log(A[3])) * exp(-(log(x / A[2]))^2 / (2log(A[3])^2))
    logn(A, x) = mapreduce((A) -> md(A, x), +, A)
    Dp = sqrt.(De[2:end] .* De[1:end-1])
    ΔlnD = log.(De[2:end] ./ De[1:end-1])
    S = logn(A, Dp)
    N = S .* ΔlnD
    return SizeDistribution(A, De, Dp, ΔlnD, S, N, :lognormal)
end

function tfun(Dl, Du, Dp)
    b = 1.0 / (-log(Dl) + log(Du))
    a = -b * log(Dl)
    return a + b * log(Dp)
end

gette(Dp) = @match Dp begin
    if Dp < 34.0
    end => 0.0
    if Dp > 1175.0
    end => 0.0
    if Dp <= 1175.0 && Dp >= 482.0
    end => tfun(1175.0, 482.0, Dp)
    if Dp <= 74.0 && Dp >= 34.0
    end => tfun(34.0, 74.0, Dp)
    _ => 1.0
end

Dp = @chain range(log10(30.0), stop = log10(1000.0), length = 1000) map(exp10, _)
𝕟 = logndist([[3224.0, 26.0, 1.97], [2168.0, 56.1, 1.39], [1262.0, 146.0, 1.44]], Dp)
df = DataFrame(Dp = Dp, t = map(gette, Dp))
ρ = 1.77
f = map(gette, 𝕟.Dp)
df1 = DataFrame(
    Dp = 𝕟.Dp .* ρ,
    S = 𝕟.S * π .* (𝕟.Dp .* 1e-3) .^ 3.0 ./ 6.0 * ρ,
    Distribution = "True",
)

df2 = DataFrame(
    Dp = 𝕟.Dp .* ρ,
    S = 𝕟.S * π .* (𝕟.Dp .* 1e-3) .^ 3.0 ./ 6.0 * ρ .* f,
    Distribution = "Atten.",
)
colors = ["black", "darkred", "steelblue3", "darkgoldenrod", "grey"]

xlabels = log10.([30, 100, 300, 1000])
xticks =
    log10.([30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
p0 = plot(
    x = 𝕟.Dp,
    y = 𝕟.S,
    Geom.line,
    Theme(default_color = "black", plot_padding = [0pt, 15pt, 0pt, 0pt]),
    Guide.xlabel("Dᵥₑ (nm)"),
    Guide.ylabel("dN/dlnD (cm⁻³)"),
    Guide.xticks(ticks = xticks),
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%i", exp10(x)) : ""),
    Scale.y_continuous(labels = x -> formatted(x, :SI, ndigits = 1)),
)

p1 = plot(
    df,
    x = :Dp,
    y = :t,
    Geom.line,
    Theme(default_color = "black", plot_padding = [-5pt, 30pt, 0pt, 0pt]),
    Guide.xlabel("Dᵥₐ (nm)"),
    Guide.ylabel("AMS fract. transmiss. (-)", orientation = :vertical),
    Guide.xticks(ticks = xticks),
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%i", exp10(x)) : ""),
)

p2 = plot(
    [df1; df2],
    x = :Dp,
    y = :S,
    color = :Distribution,
    Geom.line,
    Theme(plot_padding = [-25pt, 10pt, 0pt, 0pt]),
    Guide.xlabel("Dᵥₐ"),
    Guide.ylabel("dM/dlnD (μg m⁻³)"),
    Guide.xticks(ticks = xticks),
    Guide.colorkey(title = ""),
    Scale.color_discrete_manual(colors...),
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%i", exp10(x)) : ""),
    Coord.cartesian(xmax = 3),
)

p = hstack(p0, p1, p2)
set_default_plot_size(8inch, 2.5inch)
draw(PNG("f01.png", dpi = 600), p)
