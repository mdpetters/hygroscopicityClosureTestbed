using CSV
using DataFrames
using Dates
using Lazy
using Gadfly
using Cairo
using Fontconfig
using Measurements
using Printf
using Colors

include("TestBed.jl")

format2datetime(x) = @. DateTime(x, dateformat"mm/dd/yyyy HH:MM:SS")

function getdf(data, k1, k2)
    fields = [
        :DateTime,
        :HtoC_AMS,
        :OtoC_AMS,
        :Density_AMS,
        :Nitrate_AMS,
        :Ammonium_AMS,
        :Chloride_AMS,
        :Sulfate_AMS,
        :TotalMass_AMS,
        :OA_AMS,
        :OAFraction_AMS,
        :OADensity_ams,
    ]

    df = @> data begin
        transform(:ETime_AMS => format2datetime => :DateTime)
        transform(k1 => (x -> x) => :κ)
        transform(k2 => (x -> x) => :σκ)
        transform([:OA_AMS, :TotalMass_AMS] => ((o, t) -> o ./ t) => :OA_calc)
        transform(
            [:OA_AMS, :TotalMass_AMS] =>
                ((o, t) -> 3.0 * o ./ (2.0 .* o .+ t)) => :OA_high,
        )
        transform(
            [:OA_AMS, :TotalMass_AMS] =>
                ((o, t) -> 0.333 * o ./ (t .- 0.666 .* o)) => :OA_low,
        )
        select([fields; :OA_calc; :OA_high; :OA_low; :κ; :σκ])
        dropmissing
    end

    return df
end

function setup_case(df, i)
    n = 50
    CV = 0.05
    NH4 = df[i, :Ammonium_AMS]
    NO3 = df[i, :Nitrate_AMS]
    nrCl = df[i, :Chloride_AMS]
    SO4 = df[i, :Sulfate_AMS]
    ss = 0.0
    OA = df[i, :OA_AMS]
    rBC = 0.0
    myaw = 0.99
    tkmOA = 0.05
    tσkmOA = 0.025
    kvext = df[i, :κ]
    σkvext = df[i, :σκ]
    tρOA = df[i, :OADensity_ams]
    tσρOA = 0.2 * df[i, :OADensity_ams]

    NH4M = measurement(NH4, CV * NH4)
    SO4M = measurement(SO4, CV * SO4)
    NO3M = measurement(NO3, CV * NO3)
    nrClM = measurement(nrCl, CV * nrCl)
    ssM = measurement(ss, CV * ss)
    OAM = measurement(OA, 0.1 * OA)
    rBCM = measurement(rBC, CV * rBC)
    kmOAM = measurement(tkmOA, tσkmOA)
    ρOAM = measurement(tρOA, tσρOA)
    kvextM = measurement(kvext, σkvext)

    return TestBed.Composition(
        NH4M,
        SO4M,
        NO3M,
        nrClM,
        ssM,
        OAM,
        rBCM,
        ρOAM,
        kmOAM,
        kvextM,
        myaw,
        n,
    )
end

data = CSV.read("../data/GoAmazon2_Data_for_Markus.csv", DataFrame)
df = getdf(data, :RKappa_median_222, :Rkappa_sd_222)
λ = i -> @> setup_case(df, i) TestBed.plot_constraints
r = map(λ, 1:length(df[!, 1]))
include("plotter.jl")
