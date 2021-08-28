using DifferentialMobilityAnalyzers
using Plots
using Printf
import Plots.plot, Plots.plot!

struct OA
	ϵ::Float64
	OC::Float64
end

struct AerosolMode
	NH4::Float64
	NO3::Float64
	SO4::Float64
	nrCl::Float64
	rBC::Float64
	OA::Vector{OA}
	𝕟::SizeDistribution
	AerosolMode(NH4, NO3, SO4, nrCL, rBC, OA, 𝕟) = 
		NH4 + NO3 + SO4 + nrCL + rBC + mapfoldr(x->x.ϵ, +, OA) ≈ 1.0 ? 
		new(NH4, NO3, SO4, nrCL, rBC, OA, 𝕟) : error("Need to sum to one")
end

𝕟₁ = lognormal([[100.0, 40.0, 1.5]], d1 = 10.0, d2 = 2000.0, bins = 1024)
𝕟₂ = lognormal([[300.0, 70.0, 1.3]], d1 = 10.0, d2 = 2000.0, bins = 1024)
a = AerosolMode(0.2,0.0,0.1,0.0,0.0,[OA(0.7, 0.2)], 𝕟₁)
b = AerosolMode(0.5,0.5,0.0,0.0,0.0,[OA(0.0, 0.2)], 𝕟₂)

function getlabel(mode::AerosolMode)
	l = String[]
	mode.NH4 == 0.0 ? Nothing : push!(l, @sprintf("NH4=%0.2f ", mode.NH4))
	mode.NO3 == 0.0 ? Nothing : push!(l, @sprintf("NO3=%0.2f ", mode.NO3))
	mode.SO4 == 0.0 ? Nothing : push!(l, @sprintf("SO4=%0.2f ", mode.SO4))
	mode.nrCl == 0.0 ? Nothing : push!(l, @sprintf("SS=%0.2f ", mode.nrCl))
	mode.rBC == 0.0 ? Nothing : push!(l, @sprintf("rBC=%0.2f ", mode.rBC))
	map(1:length(mode.OA)) do i
		mode.OA[i].ϵ == 0.0 ? Nothing : push!(l, @sprintf("OA%i=%.2f ", i, mode.OA[i].ϵ))
    end
	return reduce(*, l)
end

function plot(mode::AerosolMode; kwargs...)
	plot(mode.𝕟.Dp, mode.𝕟.S; label = getlabel(mode), kwargs...)
end

function plot!(mode::AerosolMode; kwargs...)
	plot!(mode.𝕟.Dp, mode.𝕟.S; label = getlabel(mode), kwargs...)
end



plot(a; xlabel = "Volume Equivalent Diameter (nm)", ylabel = "dN/dlnD (cm-3)", 
	 xaxis = :log10, color = :red, minorgrid = true)
plot!(b, color = :green)

