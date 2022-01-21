using Pkg

Pkg.add(["Pluto", "PlutoUI", "Colors"])
Pkg.precompile()

cd("/home/pluto/notebooks")
Pkg.activate(".")
Pkg.instantiate()
