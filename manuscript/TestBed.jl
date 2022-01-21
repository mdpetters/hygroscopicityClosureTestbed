module TestBed

using Measurements
using Gadfly
using Distributions
using Random
using FileIO
using DataStructures
using DataFrames
using LinearAlgebra
using Lazy
using Printf
using Interpolations
using Statistics
using Colors
using AIOMFAC

struct Composition
    NH4::Measurement
    SO4::Measurement
    NO3::Measurement
    nrCl::Measurement
    ss::Measurement
    OA::Measurement
    rBC::Measurement
    ρOA::Measurement
    kmOA::Measurement
    kvext::Measurement
    aw::Float64
    n::Int
end

struct TestBedSolution
    ϵ::Measurement
    κm::Measurement
    κv::Measurement
    ρd::Measurement
    p1::Plot
    p2::Plot
end

function define_compounds()
    c1 = (("NH4", 2), ("SO4", 1))
    c2 = (("NH4", 1), ("NO3", 1))
    c3 = (("NH4", 1), ("Cl", 1))
    c4 = (("H", 2), ("SO4", 1))
    c5 = (("H", 1), ("NO3", 1))
    c6 = (("H", 1), ("Cl", 1))

    MWs = [132.14, 80.043, 53.491, 98.079, 63.01, 36.458]
    return Lazy.list(c1, c2, c3, c4, c5, c6), MWs
end

function define_compounds_AIOMFAC()
    c1 = (("NH4+", 2), ("SO4--", 1))
    c2 = (("NH4+", 1), ("NO3-", 1))
    c3 = (("NH4+", 1), ("Cl-", 1))
    c4 = (("H+", 2), ("SO4--", 1))
    c5 = (("H+", 1), ("NO3-", 1))
    c6 = (("H+", 1), ("Cl-", 1))

    return Lazy.list(c1, c2, c3, c4, c5, c6)
end

function create_pool(NH4, SO4, NO3, nrCl)
    nNH41 = NH4 / 18.03846
    nSO4 = SO4 / 96.0626
    nNO3 = NO3 / 62.0049
    nnrCl = nrCl / 35.453
    nH1 = 2 * nSO4 + nNO3 + nnrCl - nNH41
    nH = (nH1 < 0.0) ? 0 : nH1
    nNH4 = (nH <= 0) ? 2 * nSO4 + nNO3 + nnrCl : nNH41

    return DataFrame(H = nH, NH4 = nNH4, SO4 = nSO4, NO3 = nNO3, Cl = nnrCl)
end

function assign(pool, compounds, result)
    c = Lazy.first(compounds)
    r = Lazy.tail(compounds)
    i = @> map(x -> pool[1, x[1]] ./ x[2], c) argmin
    maxavailable = pool[1, c[i][1]] ./ c[i][2]
    n = maxavailable
    newpool = deepcopy(pool)
    newpool[1, c[1][1]] = newpool[1, c[1][1]] - n * c[1][2]
    newpool[1, c[2][1]] = newpool[1, c[2][1]] - n * c[2][2]

    if length(r) == 0
        return reverse(prepend((n, c), result))
    else
        return assign(newpool, r, prepend((n, c), result))
    end
end

function plot_AIOMFAC(df)
    set_default_plot_size(4inch, 3inch)
    ylabels = [0, 0.5, 1]
    plot(
        layer(
            df,
            x = :aw,
            y = :k,
            ymin = :kmin,
            ymax = :kmax,
            Geom.line,
            Geom.ribbon,
            Theme(
                alphas = [0.2],
                lowlight_color = c -> RGBA{Float32}(c.r, c.g, c.b),
                default_color = "black",
            ),
        ),
        Guide.xlabel("Water activity (-)"),
        Guide.ylabel("κₘ"),
        Guide.xticks(ticks = 0.5:0.1:1.0),
        Guide.yticks(ticks = 0.0:0.1:1.2),
        Scale.y_continuous(labels = y -> y in ylabels ? @sprintf("%.1f", y) : ""),
        Theme(plot_padding = [0pt, 11pt, 1pt, 0pt]),
        Coord.cartesian(xmin = 0.5, xmax = 1, ymin = 0.0, ymax = 1.2),
    )
end

function evaluate(c)
    mt = c.NH4 + c.SO4 + c.NO3 + c.nrCl + c.ss + c.OA + c.rBC
    minorg = c.NH4 + c.SO4 + c.NO3 + c.nrCl + c.ss
    finorg = minorg / mt
    fOA = c.OA / mt
    frBC = c.rBC / mt
    ρd =
        (c.OA + c.SO4 + c.NO3 + c.NH4 + c.nrCl + c.ss + c.rBC) / (
            c.OA / c.ρOA +
            (c.SO4 + c.NO3 + c.NH4) / 1.75 +
            c.nrCl / 1.52 +
            c.ss / 1.45 +
            c.rBC / 1.77
        )

    return finorg, fOA, frBC, ρd
end

function component2kappa(x)
    aw = x.RH ./ 100.0
    Gm = @. 1.0 / (1.0 - x.w)
    k = (@. (1.0 / aw - 1.0) * (Gm - 1.0))

    return aw, k
end

function computeAIOMFAC(f, MWs, ss)
    fx = define_compounds_AIOMFAC()
    fi = [0.005:0.001:0.009; 0.01:0.01:0.9]

    mass = [[f[i][1] * MWs[i] for i = 1:length(f)]; ss]
    components = [[fx[i] for i = 1:length(f)]; (("Na+", 1), ("Cl-", 1))]
    AIOMFAC.writeinput("input_0001.txt", 298.0, fi, components, mass)
    run(`./AIOMFAC Inputfiles/input_0001.txt`)
    phase, component = AIOMFAC.load_data("AIOMFAC_output_0001.txt", 298.0, fi; c = 1)
    aw, k = component2kappa(component[1])

    return DataFrame(aw = aw, k = k, comp = "Mixture")
end

function composition(x::Measurement)
    r = rand(Normal(x.val, x.err))
    out = r < 0 ? 0.0 : r

    return out
end

function monte_carlo(c)
    compounds, MWs = define_compounds()

    pool = create_pool(
        composition(c.NH4),
        composition(c.SO4),
        composition(c.NO3),
        composition(c.nrCl),
    )
    ff = assign(pool, compounds, Lazy.list())
    df = @> computeAIOMFAC(ff, MWs, c.ss.val) sort(:aw) unique
    itp = interpolate((df[!, :aw],), df[!, :k], Gridded(Linear()))

    return extrapolate(itp, 0)
end

function AIOMFAC_uncertainty(c)
    aws = @> 0.5:0.01:0.99 collect

    mat = mapfoldl(hcat, 1:c.n) do _
        extp = monte_carlo(c)
        extp(aws)
    end

    μ = mean(mat, dims = 2)[:]
    σ = std(mat, dims = 2)[:]
    i = argmin((aws .- c.aw) .^ 2.0)
    km = measurement(μ[i], σ[i])
    df = DataFrame(aw = aws, k = μ, kmin = μ .- σ, kmax = μ .+ σ)

    return df, km
end

function plot_constraints(c)
    compounds, MWs = define_compounds()

    df1, km1 = TestBed.AIOMFAC_uncertainty(c)

    OAs = 0:0.1:300
    kv1 = map(OAs) do x
        OA = measurement(x, c.OA.err)
        c1 = Composition(
            c.NH4,
            c.SO4,
            c.NO3,
            c.nrCl,
            c.ss,
            OA,
            c.rBC,
            c.ρOA,
            c.kmOA,
            c.kvext,
            c.aw,
            c.n,
        )
        fi1, fOA1, frBC1, ρd1 = evaluate(c1)
        kmo = fi1 * km1 + fOA1 * c.kmOA
        kmo * ρd1
    end

    forganic = map(OAs) do x
        OA = measurement(x, c.OA.err)
        c1 = Composition(
            c.NH4,
            c.SO4,
            c.NO3,
            c.nrCl,
            c.ss,
            OA,
            c.rBC,
            c.ρOA,
            c.kmOA,
            c.kvext,
            c.aw,
            c.n,
        )
        fi1, fOA1, frBC1, ρd1 = evaluate(c1)
        fOA1
    end

    vals = map(x -> x.val, kv1)
    err = map(x -> x.err, kv1)
    xvals = map(x -> x.val, forganic)

    i = argmin((vals .- c.kvext) .^ 2)
    mitp1 = interpolate((xvals,), vals .+ err, Gridded(Linear()))
    lvals1 = mitp1(xvals)
    i1 = argmin((lvals1 .- c.kvext.val .+ c.kvext.err) .^ 2)
    mitp2 = interpolate((xvals,), vals .- err, Gridded(Linear()))
    lvals2 = mitp2(xvals)
    i2 = argmin((lvals2 .- c.kvext.val .- c.kvext.err) .^ 2)

    p1 = plot(
        layer(
            x = xvals,
            y = vals,
            ymin = vals .- err,
            ymax = vals .+ err,
            Theme(
                alphas = [0.2],
                lowlight_color = c -> RGBA{Float32}(c.r, c.g, c.b),
                default_color = "black",
            ),
            Geom.ribbon,
            Geom.line,
        ),
        layer(
            x = [xvals[i]],
            y = [vals[i]],
            xmax = [xvals[i1]],
            xmin = [xvals[i2]],
            ymax = [c.kvext.val + c.kvext.err],
            ymin = [c.kvext.val - c.kvext.err],
            Geom.point,
            Geom.yerrorbar,
            Theme(default_color = "steelblue3"),
        ),
        layer(
            x = [xvals[i1], xvals[i2]],
            y = [c.kvext.val - c.kvext.err, c.kvext.val + c.kvext.err],
            label = ["B", "A"],
            Geom.label(position = :left),
            Geom.point,
            Theme(default_color = "steelblue3"),
        ),
        layer(
            xmax = [xvals[i1]],
            xmin = [xvals[i2]],
            ymax = [c.kvext.val + c.kvext.err],
            ymin = [c.kvext.val - c.kvext.err],
            Geom.rect,
            Theme(default_color = "steelblue3", alphas = [0.2]),
        ),
        Guide.xlabel("Organic Fraction (-)"),
        Guide.ylabel("κᵥ (-)"),
        Guide.yticks(ticks = 0:0.1:maximum(vals .+ err)),
        Guide.xticks(ticks = 0:0.1:1),
        Theme(plot_padding = [0pt, 11pt, 1pt, 0pt]),
    )

    finorg, fOA, frBC, ρd = evaluate(c)
    km = finorg * km1 + fOA * c.kmOA
    kv = km * ρd
    p2 = TestBed.plot_AIOMFAC(df1)

    ϵ = (xvals[i1] + xvals[i2]) ./ 2
    σϵ = ϵ - xvals[i2]

    return TestBedSolution(measurement(ϵ, σϵ), km, kv, ρd, p1, p2)
end

end
