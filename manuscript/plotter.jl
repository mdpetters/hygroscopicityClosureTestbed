using Statistics

ϵ = map(x -> x.κv.val, r)
σ = map(x -> x.κv.err, r)

df1 = DataFrame(t = df[!, :DateTime], kv = ϵ, kmin = ϵ .- σ, kmax = ϵ .+ σ, Color = "AMS")

df2 = DataFrame(
    t = df[!, :DateTime],
    kv = df[!, :κ],
    kmin = df[!, :κ] .- df[!, :σκ],
    kmax = df[!, :κ] .+ df[!, :σκ],
    Color = "CCN",
)

dfx = [df1; df2]

theme1 = Theme(
    alphas = [0.2],
    default_color = "steelblue3",
    line_width = 0.5pt,
    lowlight_color = c -> RGBA{Float32}(c.r, c.g, c.b),
)

colors = ["steelblue3", "darkred"]

xt = DateTime(2014, 08, 16, 0, 0, 0):Day(7):DateTime(2014, 10, 11, 0, 0, 0)
lab = DateTime(2014, 08, 16, 0, 0, 0):Day(14):DateTime(2014, 10, 11, 0, 0, 0)
p1 = plot(
    dfx,
    x = :t,
    y = :kv,
    ymin = :kmin,
    ymax = :kmax,
    color = :Color,
    Geom.line,
    Geom.ribbon,
    Guide.xlabel("Month/Day in 2014"),
    Guide.ylabel("κᵥ (-)"),
    Guide.yticks(ticks = 0:0.05:0.3),
    Guide.colorkey(""),
    Scale.x_continuous(labels = x -> x in lab ? Dates.format(x, "mm/dd") : ""),
    Scale.y_continuous(labels = y -> y in 0:0.05:0.3 ? "$y" : ""),
    Guide.xticks(ticks = xt),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(ymin = 0.0, ymax = 0.3),
    theme1,
)

datetime2time(t) = Hour(t).value #Time(Hour(t).value,  Minute(t).value, Second(t).value)

function dieldff(dfx)
    return @as df dfx begin
        transform(df, :t => (x -> datetime2time.(x)) => :Time)
        transform(df, :Time => (x -> x .+ 1) => :Time1)
        DataFrames.groupby(df, :Time)
        combine(
            df,
            :Time1 => mean => :Time1,
            :kv => mean => :kv,
            :kmin => mean => :kmin,
            :kmax => mean => :kmax,
            :Color => (x -> x[1]) => :Color,
        )
        vcat(
            df,
            DataFrame(
                Time = 24,
                Time1 = 25,
                kv = df[end, :kv],
                kmin = 0,
                kmax = 0,
                Color = df[1, :Color],
            ),
        )
    end
end

diel = [dieldff(df1); dieldff(df2)]

p3 = plot(
    diel,
    x = :Time,
    xmin = :Time,
    xmax = :Time1,
    y = :kv,
    ymin = :kmin,
    ymax = :kmax,
    color = :Color,
    Geom.rect,
    Geom.step,
    Guide.xlabel("Hour of day (UTC)"),
    Guide.ylabel("κᵥ (-)"),
    Guide.xticks(ticks = 0:1:24),
    Guide.colorkey(""),
    Guide.yticks(ticks = 0:0.05:0.3),
    Scale.x_continuous(labels = x -> x in 0:3:24 ? "$x" : ""),
    Scale.y_continuous(labels = y -> y in 0:0.05:0.3 ? "$y" : ""),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(ymin = 0, ymax = 0.3, xmax = 24),
    theme1,
)

draw(PNG("f04.png", 7inch, 3inch, dpi = 300), hstack(p1, p3))

colors = ["steelblue3", "darkred"]
ϵ = map(x -> x.ϵ.val, r)
σ = map(x -> x.ϵ.err, r)
df1 = DataFrame(t = df[!, :DateTime], kv = ϵ, kmin = ϵ .- σ, kmax = ϵ .+ σ, Color = "CCN")

df2 = DataFrame(
    t = df[!, :DateTime],
    kv = df[!, :OA_calc],
    kmin = df[!, :OA_low],
    kmax = df[!, :OA_high],
    Color = "AMS",
)

dfx = [df2; df1]
p2 = plot(
    dfx,
    x = :t,
    y = :kv,
    ymin = :kmin,
    ymax = :kmax,
    color = :Color,
    Geom.line,
    Geom.ribbon,
    Guide.ylabel("ϵₒₐ (-)"),
    Guide.yticks(ticks = 0.5:0.1:1),
    Guide.xlabel("Month/Day in 2014"),
    Guide.colorkey(""),
    Guide.xticks(ticks = xt),
    Scale.x_continuous(labels = x -> x in lab ? Dates.format(x, "mm/dd") : ""),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(ymin = 0.5, ymax = 1),
    theme1,
)

diel = [dieldff(df2); dieldff(df1)]

p4 = plot(
    diel,
    x = :Time,
    xmin = :Time,
    xmax = :Time1,
    y = :kv,
    ymin = :kmin,
    ymax = :kmax,
    color = :Color,
    Geom.step,
    Geom.rect,
    Guide.xlabel("Hour of day (UTC)"),
    Guide.ylabel("ϵₒₐ (-)"),
    Guide.colorkey(""),
    Guide.xticks(ticks = 0:1:24),
    Scale.x_continuous(labels = x -> x in 0:3:24 ? "$x" : ""),
    Guide.yticks(ticks = 0.5:0.1:1),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(ymin = 0.5, ymax = 1),
    theme1,
)

draw(PNG("f05.png", 7inch, 3inch, dpi = 300), hstack(p2, p4))
