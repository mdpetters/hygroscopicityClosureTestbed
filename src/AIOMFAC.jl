module AIOMFAC

export writeinput, parseoutput

using FileIO
using DataStructures

const header = """
Input file for AIOMFAC-web model
 
mixture components:
----
component no.:	01
component name:	'Water'
subgroup no., qty:	016, 01
----
"""

const units = """++++
mixture composition and temperature:
mass fraction?	1
mole fraction?  0	
----
"""

const group = OrderedDict([
	("Na", 202), 
	("NH4", 204), 
	("H", 205), 
	("Cl", 242), 
	("NO3", 245), 
	("HSO4", 248), 
	("SO4", 261)
])

function writeinput(file, Tk, components, fractions)
	fi = 0.001:0.01:0.9 
	open("Inputfiles/"*file, "w") do f
		write(f, header)
		map(1:length(components)) do i
			compound = components[i]
			write(f, "component no.:	$(i+1)\n")
			write(f, "component name:	'inorganic component'\n")
			map(compound) do x
				no = group[x[1]]
				write(f, "subgroup no., qty: $(no), $(x[2])\n")
			end
			write(f, "----\n")
		end
		write(f, units)
		write(f, "point, Tk, cp0, ...\n")
		map(fi) do x
			write(f, "1	    $(Tk),     ")
			map(fractions) do y
				write(f, "$(x*y/sum(fractions)),    ")
			end
			write(f,"\n")
		end
		write(f, "====")
	end

end

function mparse(x)
	return try 
		parse(Float64, x)
	catch
		0.0
	end
end

function parseoutput(file)
	lines = readlines("Outputfiles/"*file) 
	out = mapfoldl(hcat, lines[141:228]) do x
		map(mparse, split(x))
	end
	aw = out[3,:]./100.0
	w = out[4,:]
	Gm = @. 1.0/(1.0 - w)
	k = (@. (1.0/aw - 1.0) * (Gm - 1.0))
	return aw, k 
end

end
