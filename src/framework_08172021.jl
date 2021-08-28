### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ dc1a7d9c-fad9-11eb-3360-4be773edd179
begin
	using Measurements
	using Gadfly
	using PlutoUI
	using Distributions
	using Random
	using FileIO
	using DataStructures
	using DataFrames
	using ChemEquations
	using LinearAlgebra
	using Lazy
	using Printf
	using Interpolations
	using Statistics
	using Colors
	include("AIOMFAC.jl")
end

# ╔═╡ f4167b04-72b3-494e-88b5-7a7647693e2b
begin
		an = @bind n NumberField(10:10:200, default = 50)
		aCV = @bind CV NumberField(0.0:0.1:0.4, default = 0.1)
		aNH4 = @bind NH4 NumberField(0.0:0.1:2.0, default = 0.1)
		
		aNO3 = @bind NO3 NumberField(0.0:0.1:2.0, default = 0.0)
		anrCl = @bind nrCl NumberField(0.0:0.1:2.0, default = 0.0)
		aSO4 = @bind SO4 NumberField(0.0:0.1:2.0, default = 0.8)
		
		ass = @bind ss NumberField(0.0:0.1:2.0, default = 0.0)
		aOA = @bind OA NumberField(0.0:0.1:2.0, default = 0.4)
		arBC = @bind rBC NumberField(0.0:0.1:2.0, default = 0.05)
		aaw = @bind myaw NumberField(0.6:0.01:0.97, default = 0.9)
		akmOA = @bind tkmOA NumberField(0:0.01:0.2, default = 0.05)
		aσkmOA = @bind tσkmOA NumberField(0.0:0.005:0.1, default = 0.025)
		akvext = @bind kvext NumberField(0:0.01:0.7, default = 0.3)
		aσkvext = @bind σkvext NumberField(0.0:0.01:0.3, default = 0.06)
		aρOA = @bind tρOA NumberField(1.1:0.01:1.8, default = 1.5)
		aσρOA = @bind tσρOA NumberField(0.0:0.01:0.3, default = 0.1)
	
		kmOA = measurement(tkmOA, tσkmOA) 
		ρOA = measurement(tρOA, tσρOA) 
		
	md"""**Use the number fields to set the values**
		
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

# ╔═╡ a522db23-4c9a-4964-90eb-fc6c6a54da18
begin
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
	
	compounds, MWs = define_compounds()
end

# ╔═╡ 08c6d14d-7dcd-4025-bf28-186a90bd7cd5
function create_pool(NH4, SO4, NO3, nrCl, ss)	
	nNH41 = NH4/18.03846
	nSO4 = SO4/96.0626
	nNO3 = NO3/62.0049
	nnrCl = nrCl/35.453
	nH1 = 2*nSO4 + nNO3 + nnrCl - nNH41
	nH = (nH1 < 0.0) ? 0 : nH1
	nNH4 = (nH <= 0) ? 2*nSO4 + nNO3 + nnrCl : nNH41
		
	return DataFrame(
		H = nH, 
		NH4 = nNH4, 
		Na = ss/22.989769, 
		SO4 = nSO4, 
		NO3 = nNO3,
		Cl = nnrCl + ss/35.453
	)
end


# ╔═╡ f6cab4cb-cf15-4d2f-9056-b960343ee83f
begin
	mypool = create_pool(NH4, SO4, NO3, nrCl, ss)

md"""**Molar composition**
	
	Computed molar ratios in units of μmol. If the mass fractions are selected such that there is an excess of cations, implying that the H+ concentration is < 0, then the NH4+ concentration is set to balance the anions. The shown molar composition is used to initialize the AIOMFAC model.
	
	$(mypool)	
	

	"""
end


# ╔═╡ 992ba913-a7fd-4bfb-a7e6-b943fd422421
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

# ╔═╡ 53048c16-d8c5-40b5-b824-0e6ccf598678
function plot_AIOMFAC(df)
	set_default_plot_size(4inch, 3inch)
	ylabels = [0, 0.5, 1]
	plot(layer(df, x=:aw, y=:k, ymin = :kmin, ymax = :kmax, Geom.line, Geom.ribbon,
		Theme(
			alphas = [0.2], 
			lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b), 
			default_color="black"
		)),
		Guide.xlabel("Water activity (-)"),
		Guide.ylabel("κₘ"),
		Guide.xticks(ticks=0.5:0.1:1.0),
		Guide.yticks(ticks=0.0:0.1:1.2),
	    Scale.y_continuous(labels=y -> y in ylabels ? @sprintf("%.1f", y) : ""),
		Theme(plot_padding=[0pt, 11pt, 1pt, 0pt]),
		Coord.cartesian(xmin=0.5, xmax=1, ymin=0.0, ymax=1.2)
	)
end


# ╔═╡ d394eaf9-2f4f-4f43-893e-1bb05ac905dc
function evaluate(NH4, SO4, NO3, nrCl, ss, OA, rBC, ρOA)
	NH4M = measurement(NH4, CV*NH4)
	SO4M = measurement(SO4, CV*SO4)
	NO3M = measurement(NO3, CV*NO3)
	nrClM = measurement(nrCl, CV*nrCl)
	ssM = measurement(ss, CV*ss)
	OAM = measurement(OA, CV*OA)
	rBCM = measurement(rBC, CV*rBC)
	mt = NH4M + SO4M + NO3M + nrClM + ssM + OAM + rBCM 
	minorg = NH4M + SO4M + NO3M + nrClM + ssM 
	finorg = minorg/mt
	fOA = OA/mt
	frBC = rBC/mt
	ρd = (OAM + SO4M + NO3M + NH4M + nrClM + ssM + rBCM)/
		(OAM/ρOA + (SO4M + NO3M + NH4M)/1.75 + nrClM/1.52 + ssM/1.45 + rBCM/1.77)
	return finorg, fOA, frBC, ρd
end

# ╔═╡ c98fc5de-a111-4d91-b7fa-056bb1847093
function computeAIOMFAC(f, MWs, ss)	
	mass = [[f[i][1]*MWs[i] for i = 1:length(f)]; ss]
	components = [[f[i][2] for i = 1:length(f)]; (("Na", 1), ("Cl", 1))]
	AIOMFAC.writeinput("input_0001.txt", 298, components, mass)
	run(`./AIOMFAC Inputfiles/input_0001.txt`)
	aw, k = AIOMFAC.parseoutput("AIOMFAC_output_0001.txt")
	return DataFrame(aw=aw, k=k, comp="Mixture")
end

# ╔═╡ 2bedbea1-3aee-4b6c-9c4c-1f68b66b7a76
function composition(μ, σ)
	r = rand(Normal(μ, σ))
	out = r < 0 ? 0.0 : r
	return out
end

# ╔═╡ 0a3a331a-805b-45dc-bf92-0a59cf89e3ce
function monte_carlo(compounds, MWs)
	pool = create_pool(
		composition(NH4, NH4*CV), 
		composition(SO4, SO4*CV), 
		composition(NO3, NO3*CV), 
		composition(nrCl, nrCl*CV),
		composition(ss, ss*CV)
	)
	ff = assign(pool, compounds, Lazy.list())
	df = computeAIOMFAC(ff, MWs, ss)
	itp = interpolate((reverse(df[!,:aw]),), reverse(df[!,:k]), Gridded(Linear()))
	extrapolate(itp, 0)
end

# ╔═╡ bc04839f-9807-485f-b1e6-256ce8f9deff
function AIOMFAC_uncertainty(compounds, MWs, myaw, n)
	aws = @> 0.5:0.01:0.98 collect
	mat = mapfoldl(hcat, 1:n) do _
		extp = monte_carlo(compounds, MWs) 
		extp(aws)
	end
	μ = mean(mat, dims = 2)[:]
	σ = std(mat, dims = 2)[:]
	i = argmin((aws .- myaw).^2.0)
	km = measurement(μ[i], σ[i])
	df = DataFrame(aw = aws, k = μ, kmin = μ.-σ, kmax = μ.+σ)
	return df, km
end

# ╔═╡ d6f4dabe-1f10-4a2a-95ff-fafd1bb17127
begin
	
	fx = assign(mypool, compounds, Lazy.list())
	df1, km1 = AIOMFAC_uncertainty(compounds, MWs, myaw, n)
	p = plot_AIOMFAC(df1)

md"""**AIMOFAC Output**
	
	Computed κₘ vs. water activity from AIOMFAC simulation for molar composition displayed above.
	
	$p

	"""
end


# ╔═╡ 06919bc9-6135-4c63-a7cb-559833501af9
begin
	finorg, fOA, frBC, ρd = evaluate(NH4, SO4, NO3, nrCl, ss, OA, rBC, ρOA)
	md"""
	**Evaluated Uncertainties**
	
	The error for finorg, forg is determined using error propagation. The error in km,org is specified. The error in km,inorg is determined via Monte Carlo simulations with randomly perturbed inputs to AIOMFAC. The error in density is computed using error propagation.
	
	 | f ± σ | κₘ ± σ | ρd ± σ |
	| ----------- | ----------- |----------- |----------- |
	|Inorganic | $(finorg) | $(km1) | $ρd |
	|Organic  | $(fOA) | $(kmOA) | |
	| rBC  | $(frBC) | 0 ± 0 | |
	
	"""
end

# ╔═╡ c6eba237-da8c-4b66-add5-39a96df897ad
begin 
	km = finorg*km1 + fOA*kmOA
	kv = km*ρd

	md"""
	**AMS-derived Based Hygroscopcity Parameters**
	
	The value is computed from the table above using error propagation.
	
	$$\kappa_m = \sum f_i \kappa_{m,i}$$ 
	
	||
	|:--:| 
	| κₘ = $(km)  | 

	$$\kappa_v = \kappa_{m} \frac{ρ_d}{ρ_w}$$ 
	
	||
	|:--:| 
	| κᵥ = $(kv)  | 
	
	"""
end

# ╔═╡ 5db9ae42-8c57-490e-a4ee-eb650f78c9cb
begin
	OAs = 0:0.1:300
	kv1 = map(OAs) do x
		fi1, fOA1, frBC1, ρd1 = evaluate(NH4, SO4, NO3, nrCl, ss, x, rBC, ρOA)
		kmo = fi1*km1 + fOA1*kmOA 
		kmo*ρd1
	end
	
	forganic = map(OAs) do x
		fi1, fOA1, frBC1, ρd1 = evaluate(NH4, SO4, NO3, nrCl, ss, x, rBC, ρOA)
		fOA1
	end
		
	vals = map(x->x.val, kv1)
	err = map(x->x.err, kv1)
	xvals = map(x->x.val, forganic)
	
	i = argmin((vals .- kvext).^2)
	mitp1 = interpolate((xvals,),vals .+ err, Gridded(Linear()))
	lvals1 = mitp1(xvals)
	i1 = argmin((lvals1 .- kvext .+ σkvext).^2)
	mitp2 = interpolate((xvals,),vals .- err, Gridded(Linear()))
	lvals2 = mitp2(xvals)
	i2 = argmin((lvals2 .- kvext .- σkvext).^2)
	
	set_default_plot_size(4inch, 3inch)
	p2 = plot(
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
			ymax = [kvext + σkvext],
			ymin = [kvext - σkvext],
			Geom.point, 
			Geom.yerrorbar, 
			Theme(default_color = "steelblue3")
		),
		layer(
			x = [xvals[i1], xvals[i2]],
 			y = [kvext - σkvext, kvext + σkvext],
			label = ["B", "A"],
			Geom.label(position=:left),
			Geom.point,
			Theme(default_color = "steelblue3")
		),
		layer(
			xmax = [xvals[i1]], 
			xmin = [xvals[i2]],
			ymax = [kvext + σkvext],
			ymin = [kvext - σkvext],
			Geom.rect, 
			Theme(default_color = "steelblue3", alphas = [0.2])
		),
		Guide.xlabel("Organic Fraction (-)"),
		Guide.ylabel("κᵥ (-)"),
		Guide.yticks(ticks = 0:0.1:maximum(vals .+ err)),
		Guide.xticks(ticks = 0:0.1:1)
	)
	
	
	md"""
	**External Constraint on OA mass**
	
	Assuming that κᵥ±σₖᵥ is measured using a second instrument, what are the valid ranges of organic aerosol mass (or AMS-derived organic fraction) that will provide closure the AMS-derived κᵥ?
	
	$p2	
	
	If the externally measured κᵥ = $(measurement(kvext,σkvext)) then  $(round(xvals[i2], digits = 2)) < forg < $(round(xvals[i1], digits = 2)) are valid within measurement uncertainties. This corresponds to $(round(OAs[i2], digits = 2)) < $(OAs[i]) < $(round(OAs[i1], digits = 2)) are valid within measurement uncertainties.



	"""
	
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ChemEquations = "75e6b969-eb85-49eb-b654-b2cb515e226f"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
Lazy = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
ChemEquations = "~0.2.1"
Colors = "~0.12.8"
DataFrames = "~1.2.2"
DataStructures = "~0.18.10"
Distributions = "~0.24.18"
FileIO = "~1.10.1"
Gadfly = "~1.3.3"
Interpolations = "~0.13.4"
Lazy = "~0.15.1"
Measurements = "~2.6.0"
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "a4d07a1c313392a77042855df46c5f534076fab9"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[CategoricalArrays]]
deps = ["DataAPI", "Future", "JSON", "Missings", "Printf", "RecipesBase", "Statistics", "StructTypes", "Unicode"]
git-tree-sha1 = "1562002780515d2573a4fb0c3715e4e57481075e"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "bdc0937269321858ab2a4f288486cb258b9a0af7"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.3.0"

[[ChemEquations]]
deps = ["DocStringExtensions", "LinearAlgebraX"]
git-tree-sha1 = "9dcd2d290940c9a1f40c69414dc567b6ddf6f73a"
uuid = "75e6b969-eb85-49eb-b654-b2cb515e226f"
version = "0.2.1"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "344f143fa0ec67e47917848795ab19c6a455f32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.32.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "c6461fc7c35a4bb8d00905df7adafcff1fe3a6bc"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.2"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[CoupledFields]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "6c9671364c68c1158ac2524ac881536195b7e7bc"
uuid = "7ad07ef1-bdf2-5661-9d2b-286fd4296dac"
version = "0.2.0"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d785f42445b63fc86caa08bb9a9351008be9b765"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.2.2"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "abe4ad222b26af3337262b8afb28fab8d215e9f8"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.3"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "a837fdf80f333415b69684ba8e8ae6ba76de6aaa"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.24.18"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f985af3b9f4e278b1d24434cbb546d6092fca661"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.3"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3676abafff7e4ff07bbd2c42b3d8201f31653dcc"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.9+8"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "256d8e6188f3f1ebfa1a5d17e072a0efafa8c5bf"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.10.1"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "693210145367e7685d8604aee33d9bfb85db8b31"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.11.9"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[Gadfly]]
deps = ["Base64", "CategoricalArrays", "Colors", "Compose", "Contour", "CoupledFields", "DataAPI", "DataStructures", "Dates", "Distributions", "DocStringExtensions", "Hexagons", "IndirectArrays", "IterTools", "JSON", "Juno", "KernelDensity", "LinearAlgebra", "Loess", "Measures", "Printf", "REPL", "Random", "Requires", "Showoff", "Statistics"]
git-tree-sha1 = "96da4818e4d481a29aa7d66aac1eb778432fb89a"
uuid = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
version = "1.3.3"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[Hexagons]]
deps = ["Test"]
git-tree-sha1 = "de4a6f9e7c4710ced6838ca906f81905f7385fd6"
uuid = "a1b4810d-1bce-5fbd-ac56-80944d57a21f"
version = "0.2.0"

[[IndirectArrays]]
git-tree-sha1 = "c2a145a145dc03a7620af1444e0264ef907bd44f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "0.5.1"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "61aa005707ea2cebf47c8d780da8dc9bc4e0c512"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.4"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[InvertedIndices]]
deps = ["Test"]
git-tree-sha1 = "15732c475062348b0165684ffe28e85ea8396afc"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.0.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

[[KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LinearAlgebraX]]
deps = ["LinearAlgebra", "Mods", "SimplePolynomials"]
git-tree-sha1 = "b7d3f84cad6c72c13c4ce8d705f1fa3f49dede22"
uuid = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
version = "0.0.6"

[[Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "b5254a86cf65944c68ed938e575f5c81d5dfe4cb"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.5.3"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "7bd5f6565d80b6bf753738d2bc40a5dfea072070"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.2.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "c253236b0ed414624b083e6b72bfe891fbd2c7af"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+1"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "0fb723cd8c45858c22169b2e42269e53271a6df7"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.7"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "31c8c0569b914111c94dd31149265ed47c238c5b"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.6.0"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "4ea90bd5d3985ae1f9a908bd4500ae88921c5ce7"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "748f6e1e4de814b101911e64cc12d83a6af66782"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.2"

[[Mods]]
git-tree-sha1 = "75ef35a0a13acf77a6a1a098af6607f810d36780"
uuid = "7475f97c-0381-53b1-977b-4c60186c8d62"
version = "1.3.1"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[Multisets]]
git-tree-sha1 = "3c478ef38e8d858aed1aeba3a2043be72154e3c7"
uuid = "3b2b4ff1-bcff-5658-a3ee-dbcf1ce5ac09"
version = "0.4.2"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "c0f4a4836e5f3e0763243b8324200af6d0e0f90c"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.5"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "477bf42b4d1496b454c10cce46645bb5b8a0cf2c"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.2"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "OffsetArrays", "RecipesBase"]
git-tree-sha1 = "0b15f3597b01eb76764dd03c3c23d6679a3c32c8"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "1.2.1"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "cde4ce9d6f33219465b55162811d8de8139c0414"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.2.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "0d1245a357cc61c8cd61934c07447aa569ff22e6"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.1.0"

[[Primes]]
git-tree-sha1 = "afccf037da52fa596223e5a0e331ff752e0e845c"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "7dff99fbc740e2f8228c6878e2aad6d7c2678098"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.1"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimplePolynomials]]
deps = ["Mods", "Multisets", "Polynomials", "Primes"]
git-tree-sha1 = "12e65afc86adb18265de033311b442b2575747ad"
uuid = "cc47b68c-3164-5771-a705-2bc0097375a0"
version = "0.2.7"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "a322a9493e49c5f3a10b50df3aedaf1cdb3244b7"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "fed1ec1e65749c4d96fc20dd13bea72b55457e62"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.9"

[[StatsFuns]]
deps = ["LogExpFunctions", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "30cd8c360c54081f806b1ee14d2eecbef3c04c49"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.8"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "e36adc471280e8b346ea24c5c87ba0571204be7a"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.7.2"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "d0c690d37c73aeb5ca063056283fde5585a41710"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "81753f400872e5074768c9a77d4c44e70d409ef0"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.6"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "59e2ad8fd1591ea019a5259bd012d7aee15f995c"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─dc1a7d9c-fad9-11eb-3360-4be773edd179
# ╟─f4167b04-72b3-494e-88b5-7a7647693e2b
# ╟─f6cab4cb-cf15-4d2f-9056-b960343ee83f
# ╟─d6f4dabe-1f10-4a2a-95ff-fafd1bb17127
# ╟─06919bc9-6135-4c63-a7cb-559833501af9
# ╟─c6eba237-da8c-4b66-add5-39a96df897ad
# ╟─5db9ae42-8c57-490e-a4ee-eb650f78c9cb
# ╟─a522db23-4c9a-4964-90eb-fc6c6a54da18
# ╟─08c6d14d-7dcd-4025-bf28-186a90bd7cd5
# ╟─992ba913-a7fd-4bfb-a7e6-b943fd422421
# ╟─53048c16-d8c5-40b5-b824-0e6ccf598678
# ╟─d394eaf9-2f4f-4f43-893e-1bb05ac905dc
# ╟─c98fc5de-a111-4d91-b7fa-056bb1847093
# ╟─0a3a331a-805b-45dc-bf92-0a59cf89e3ce
# ╟─bc04839f-9807-485f-b1e6-256ce8f9deff
# ╟─2bedbea1-3aee-4b6c-9c4c-1f68b66b7a76
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
