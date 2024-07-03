# EmbraceUncertainty

The book "Embrace Uncertainty: Fitting Mixed-Effects Models with Julia"

This repository uses [Quarto](https://quarto.org) with the Julia code execution supplied by [QuartoNotebookRunner.jl](https://github.com/PumasAI/QuartoNotebookRunner.jl/), which requires Quarto 1.5+.

```sh
~/EmbraceUncertainty$ julia

julia> using Pkg

julia> Pkg.activate(".")
  Activating project at `~/EmbraceUncertainty`

julia> Pkg.instantiate()
< lots of output >

julia> exit()

~/EmbraceUncertainty$ quarto preview

< lots of output >
```
