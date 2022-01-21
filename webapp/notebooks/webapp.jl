### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ dc1a7d9c-fad9-11eb-3360-4be773edd179
begin    
	import Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

	using Measurements
	using Gadfly
	using PlutoUI
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

	md"""
	# Hygroscopicity Closure Testbed 

	Welcome to this interactive testbed application. Use the number fields to change the composition of the aerosol. After each change, the reactive Pluto notebook revaluates all of the cells. This behavior is equivalent to that of a spreadsheet application.  

	If you have questions or comments, please send an email to\
	Markus Petters\
	**mdpetter@ncsu.edu**

	**Source Code** 
	
	Pluto notebooks are referentially transparent. The complete source code of the testbed is included in the next cell. Use show/hide code icon to the left view the code. 
	"""
end

# ╔═╡ 5c8493d3-8217-45c2-8667-23848715cdb8
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
	nNH41 = NH4/18.03846
	nSO4 = SO4/96.0626
	nNO3 = NO3/62.0049
	nnrCl =nrCl/35.453
	nH1 = 2*nSO4 + nNO3 + nnrCl - nNH41
	nH = (nH1 < 0.0) ? 0 : nH1
	nNH4 = (nH <= 0) ? 2*nSO4 + nNO3 + nnrCl : nNH41

	return DataFrame(
					 H = nH, 
					 NH4 = nNH4, 
					 SO4 = nSO4, 
					 NO3 = nNO3,
					 Cl = nnrCl 
					 )
end

function assign(pool, compounds, result)
	c = Lazy.first(compounds)
	r = Lazy.tail(compounds)
	i = @> map(x->pool[1,x[1]]./x[2], c) argmin
	maxavailable = pool[1,c[i][1]]./c[i][2]
	n = maxavailable
	newpool = deepcopy(pool)
	newpool[1,c[1][1]] = newpool[1,c[1][1]] - n*c[1][2]
	newpool[1,c[2][1]] = newpool[1,c[2][1]] - n*c[2][2]
	
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
		 layer(df, x=:aw, y=:k, ymin = :kmin, ymax = :kmax, Geom.line, Geom.ribbon,
			   Theme(
					 alphas = [0.2], 
					 lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b), 
					 default_color="black"
					 )
			   ),
		 Guide.xlabel("Water activity (-)"),
		 Guide.ylabel("κₘ"),
		 Guide.xticks(ticks=0.5:0.1:1.0),
		 Guide.yticks(ticks=0.0:0.1:1.2),
		 Scale.y_continuous(labels=y -> y in ylabels ? @sprintf("%.1f", y) : ""),
		 Theme(plot_padding=[0pt, 11pt, 1pt, 0pt]),
		 Coord.cartesian(xmin=0.5, xmax=1, ymin=0.0, ymax=1.2)
		 )
end

function evaluate(c)
	mt = c.NH4 + c.SO4 + c.NO3 + c.nrCl + c.ss + c.OA + c.rBC 
	minorg = c.NH4 + c.SO4 + c.NO3 + c.nrCl + c.ss 
	finorg = minorg/mt
	fOA = c.OA/mt
	frBC = c.rBC/mt
	ρd = (c.OA + c.SO4 + c.NO3 + c.NH4 + c.nrCl + c.ss + c.rBC)/
	(c.OA/c.ρOA + (c.SO4 + c.NO3 + c.NH4)/1.75 + c.nrCl/1.52 + c.ss/1.45 + c.rBC/1.77)

	return finorg, fOA, frBC, ρd
end

function component2kappa(x)
	aw = x.RH./100.0
	Gm = @. 1.0/(1.0 - x.w)
	k = (@. (1.0/aw - 1.0) * (Gm - 1.0))

	return aw, k
end

function computeAIOMFAC(f, MWs, ss)	
	fx = define_compounds_AIOMFAC()
	fi = [0.005:0.001:0.009;0.01:0.01:0.9]

	mass = [[f[i][1]*MWs[i] for i = 1:length(f)]; ss]
	components = [[fx[i] for i = 1:length(f)]; (("Na+", 1), ("Cl-", 1))]
	AIOMFAC.writeinput("input_0001.txt", 298.0, fi, components, mass)
	run(`./AIOMFAC Inputfiles/input_0001.txt`)
	phase, component = AIOMFAC.load_data("AIOMFAC_output_0001.txt", 298.0, fi;c=1)
	aw, k = component2kappa(component[1]) 

	return DataFrame(aw=aw, k=k, comp="Mixture")
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
	knots = Interpolations.deduplicate_knots!(df[!,:aw])
	itp = interpolate((knots,), df[!,:k], Gridded(Linear()))

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
	i = argmin((aws .- c.aw).^2.0)
	km = measurement(μ[i], σ[i])
	df = DataFrame(aw = aws, k = μ, kmin = μ.-σ, kmax = μ.+σ)

	return df, km
end

function plot_constraints(c)
	compounds, MWs = define_compounds()

	df1, km1 = TestBed.AIOMFAC_uncertainty(c)
	
	OAs = 0:0.1:300
	kv1 = map(OAs) do x
		OA = measurement(x, c.OA.err)
		c1 = Composition(c.NH4, c.SO4, c.NO3, c.nrCl, c.ss, OA, c.rBC, c.ρOA, c.kmOA, 
						 c.kvext, c.aw, c.n)
		fi1, fOA1, frBC1, ρd1 = evaluate(c1)
		kmo = fi1*km1 + fOA1*c.kmOA 
		kmo*ρd1
	end

	forganic = map(OAs) do x
		OA = measurement(x, c.OA.err)
		c1 = Composition(c.NH4, c.SO4, c.NO3, c.nrCl, c.ss, OA, c.rBC, c.ρOA, c.kmOA, 
						 c.kvext, c.aw, c.n)
		fi1, fOA1, frBC1, ρd1 = evaluate(c1)
		fOA1
	end

	vals = map(x->x.val, kv1)
	err = map(x->x.err, kv1)
	xvals = map(x->x.val, forganic)

	i = argmin((vals .- c.kvext).^2)
	knots = Interpolations.deduplicate_knots!(xvals)
	mitp1 = interpolate((knots,),vals .+ err, Gridded(Linear()))
	lvals1 = mitp1(xvals)
	i1 = argmin((lvals1 .- c.kvext.val .+ c.kvext.err).^2)
	mitp2 = interpolate((knots,),vals .- err, Gridded(Linear()))
	lvals2 = mitp2(xvals)
	i2 = argmin((lvals2 .- c.kvext.val .- c.kvext.err).^2)

	p1 = plot(
			  layer(x = xvals, y = vals, ymin = vals .- err, ymax = vals .+ err, 
					Theme(
						  alphas = [0.2], 
						  lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b), 
						  default_color="black"
						  ),
					Geom.ribbon, 
					Geom.line),
			  layer(
					x = [xvals[i]], 
					y = [vals[i]], 
					xmax = [xvals[i1]], 
					xmin = [xvals[i2]],
					ymax = [c.kvext.val + c.kvext.err],
					ymin = [c.kvext.val - c.kvext.err],
					Geom.point, 
					Geom.yerrorbar, 
					Theme(default_color = "steelblue3")
					),
			  layer(
					x = [xvals[i1], xvals[i2]],
					y = [c.kvext.val - c.kvext.err, c.kvext.val + c.kvext.err],
					label = ["B", "A"],
					Geom.label(position=:left),
					Geom.point,
					Theme(default_color = "steelblue3")
					),
			  layer(
					xmax = [xvals[i1]], 
					xmin = [xvals[i2]],
					ymax = [c.kvext.val + c.kvext.err],
					ymin = [c.kvext.val - c.kvext.err],
					Geom.rect, 
					Theme(default_color = "steelblue3", alphas = [0.2])
					),
			  Guide.xlabel("Organic Fraction (-)"),
			  Guide.ylabel("κᵥ (-)"),
			  Guide.yticks(ticks = 0:0.1:maximum(vals .+ err)),
			  Guide.xticks(ticks = 0:0.1:1),
			  Theme(plot_padding=[0pt, 11pt, 1pt, 0pt]),
			  )

	finorg, fOA, frBC, ρd = evaluate(c)
	km = finorg*km1 + fOA*c.kmOA
	kv = km*ρd
	p2 = TestBed.plot_AIOMFAC(df1)

	ϵ = (xvals[i1] + xvals[i2])./2
	σϵ = ϵ - xvals[i2]

	return TestBedSolution(measurement(ϵ,σϵ), km, kv, ρd, p1, p2)
end

end

# ╔═╡ f4167b04-72b3-494e-88b5-7a7647693e2b
begin
	an = @bind n NumberField(10:10:200, default = 20)
	aCV = @bind CV NumberField(0.0:0.1:0.4, default = 0.1)
	aNH4 = @bind NH4 NumberField(0.0:0.1:2.0, default = 0.1)	
	aNO3 = @bind NO3 NumberField(0.0:0.1:2.0, default = 0.0)
	anrCl = @bind nrCl NumberField(0.0:0.1:2.0, default = 0.0)
	aSO4 = @bind SO4 NumberField(0.0:0.1:2.0, default = 0.8)
	ass = @bind ss NumberField(0.0:0.1:2.0, default = 0.0)
	aOA = @bind OA NumberField(0.0:0.1:4.0, default = 0.4)
	arBC = @bind rBC NumberField(0.0:0.1:2.0, default = 0.05)
	aaw = @bind myaw NumberField(0.6:0.01:0.97, default = 0.9)
	akmOA = @bind tkmOA NumberField(0:0.01:0.2, default = 0.05)
	aσkmOA = @bind tσkmOA NumberField(0.0:0.005:0.1, default = 0.025)
	akvext = @bind kvext NumberField(0:0.01:0.7, default = 0.3)
	aσkvext = @bind σkvext NumberField(0.0:0.01:0.3, default = 0.06)
	aρOA = @bind tρOA NumberField(1.1:0.01:1.8, default = 1.5)
	aσρOA = @bind tσρOA NumberField(0.0:0.01:0.3, default = 0.1)
	
	md"""
		All concentrations are in μg m⁻³, CV is the coefficient of variation for σᵢ/μᵢ for all mass concentrations, n is the number of simulations, aw ± σaw is target water activity and it's uncertainty, κm,OA ± σκm,OA is the hygroscopicity parameter of the organic fraction, and κv,ext ± σκv,ext is the volume-based hygroscopicity parameter from the external constraint. 
		
		| Problem Setup      | Cations (μg m⁻³) | Anions (μg m⁻³) | Other (μg m⁻³) | 
		| ----------- | ----------- | ----------- | ----------- | 
		| CV : $(aCV)    | NH₄⁺: $(aNH4)| SO₄²⁻: $(aSO4)	| Seasalt: $(ass) | 
		| n: $(an) | | NO₃⁻: $(aNO3)| OA: $(aOA)	 | 
		| aw: $(aaw) |  | nrCl⁻: $(anrCl) | rBC: $(arBC) | 
		| ρOA,ext: $(aρOA) | | |
		| σρOA,ext: $(aσρOA) | | | |
		| κm,OA: $(akmOA) | | |
		| σκm,OA: $(aσkmOA) | | | |
		| κv,ext: $(akvext) | | |
		| σκv,ext: $(aσkvext) | | | |
		
		
		"""
	
end

# ╔═╡ 9ad20d4a-21e0-491e-b7fc-9c5d10a9bc45
begin
	NH4M = measurement(NH4, CV*NH4)
	SO4M = measurement(SO4, CV*SO4)
	NO3M = measurement(NO3, CV*NO3)
	nrClM = measurement(nrCl, CV*nrCl)
	ssM = measurement(ss, CV*ss)
	OAM = measurement(OA, 0.1*OA)
	rBCM = measurement(rBC, CV*rBC)
	kmOAM = measurement(tkmOA, tσkmOA) 
	ρOAM = measurement(tρOA, tσρOA) 
	kvextM = measurement(kvext, σkvext)

	comp = TestBed.Composition(NH4M, SO4M, NO3M, nrClM, ssM, OAM, rBCM, ρOAM, 
							   kmOAM, kvextM, myaw, n)
	sol = TestBed.plot_constraints(comp)

	compounds, MWs = TestBed.define_compounds()

	mypool = TestBed.create_pool(
					   TestBed.composition(comp.NH4), 
					   TestBed.composition(comp.SO4), 
					   TestBed.composition(comp.NO3), 
					   TestBed.composition(comp.nrCl),
					   )
	ff = TestBed.assign(mypool, compounds, Lazy.list())
	a = map(x -> @sprintf("\n%f, %s", x[1], x[2]), ff);

md"""**Molar composition**
	
	Computed molar ratios in units of μmol. If the mass fractions are selected such that there is an excess of cations, implying that the H+ concentration is < 0, then the NH4+ concentration is set to balance the anions. The shown molar composition is used to initialize the AIOMFAC model.
	
	$(mypool)	

	The initialized composition passed to AIOMFAC is

	"""
end


# ╔═╡ d6f4dabe-1f10-4a2a-95ff-fafd1bb17127
begin

set_default_plot_size(7inch, 3inch) 
	
md"""**AIMOFAC Output**
	
	**Left:** Computed κₘ vs. water activity from AIOMFAC simulation for molar composition displayed above. **Right:** Modeled kv versus organic fraction. In both plots, the grey shaded area corresponds to the uncertainty derived from the Monte-Carlo simulation. The blue vertical error bar shows the external constraint kvext.  the points A and B correspond to lowest and highest organic fraction that is consistent with the external constraint. 

	$(hstack(sol.p2, sol.p1))

	Assuming that κᵥ±σₖᵥ is measured using a second instrument, what are the valid ranges of organic aerosol mass (or AMS-derived organic fraction) that will provide closure the AMS-derived κᵥ? If the externally measured κᵥ = $(kvextM) then $(round(sol.ϵ.val-sol.ϵ.err, digits = 2)) < ϵ < $(round(sol.ϵ.val+sol.ϵ.err, digits = 2)) are valid within measurement uncertainties. 
	"""

end


# ╔═╡ 85139c4e-b6a7-43da-9a5a-6443400e4223
begin
	df1, km1 = TestBed.AIOMFAC_uncertainty(comp) 
	finorg, fOA, frBC, ρd = TestBed.evaluate(comp)
	
	md"""
	**Evaluated Uncertainties**
	
	The error for finorg, forg is determined using error propagation. The error in km,org is specified. The error in km,inorg is determined via Monte Carlo simulations with randomly perturbed inputs to AIOMFAC. The error in density is computed using error propagation.
	
	 | ϵ ± σ | κₘ ± σ | ρd ± σ |
	| ----------- | ----------- |----------- |----------- |
	|Inorganic | $(finorg) | $(km1) | $ρd |
	|Organic  | $(fOA) | $(kmOAM) | |
	| rBC  | $(frBC) | 0 ± 0 | |
	
	"""
end

# ╔═╡ 5c526ad1-902c-4ee5-8b6f-d786d36e3ef3
begin 
	km = finorg*km1 + fOA*kmOAM
	kv = km*ρd

	md"""
	**AMS-derived Based Hygroscopcity Parameters**
	
	The value is computed from the table above using error propagation.
	
	$$\kappa_m = \sum \epsilon_i \kappa_{m,i}$$ 
	
	||
	|:--:| 
	| κₘ = $(km)  | 

	$$\kappa_v = \kappa_{m} \frac{ρ_d}{ρ_w}$$ 
	
	||
	|:--:| 
	| κᵥ = $(kv)  | 
	
	"""
end

# ╔═╡ 9858f59b-8c79-4f2b-8dc8-e9976366f14d
a

# ╔═╡ Cell order:
# ╟─dc1a7d9c-fad9-11eb-3360-4be773edd179
# ╟─5c8493d3-8217-45c2-8667-23848715cdb8
# ╟─f4167b04-72b3-494e-88b5-7a7647693e2b
# ╟─d6f4dabe-1f10-4a2a-95ff-fafd1bb17127
# ╟─9ad20d4a-21e0-491e-b7fc-9c5d10a9bc45
# ╟─5c526ad1-902c-4ee5-8b6f-d786d36e3ef3
# ╟─85139c4e-b6a7-43da-9a5a-6443400e4223
# ╟─9858f59b-8c79-4f2b-8dc8-e9976366f14d
