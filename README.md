# EmbraceUncertainty

The book "Embrace Uncertainty: Fitting Mixed-Effects Models with Julia"

This repository uses [Quarto](https://quarto.org). To be able to render all the pages, you will need an appropriate Jupyter kernel installed and the local environment instantiated.

```sh
~/EmbraceUncertainty$ julia

julia> using Pkg

julia> Pkg.add("IJulia")
    Updating registry at `~/.julia/registries/General.toml`
   Resolving package versions...
< lots of output >

julia> using IJulia

julia> installkernel("julia", "--threads=auto", "--project=@.")
[ Info: Installing julia kernelspec in ~/.local/share/jupyter/kernels/julia-1.9
"~/.local/share/jupyter/kernels/julia-1.9"

julia> Pkg.activate(".")
  Activating project at `~/EmbraceUncertainty`

julia> Pkg.instantiate()
< lots of output >

julia> exit()

~/EmbraceUncertainty$ quarto preview

< lots of output >

```

Note two important parts of the kernel used above: a flag to activate local Julia environment and a flag to allow Julia to automatically set the appropriate number of threads.
The former is necessary to render the book at all; the latter may improve render speed.
