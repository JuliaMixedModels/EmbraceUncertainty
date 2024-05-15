# EmbraceUncertainty

The book "Embrace Uncertainty: Fitting Mixed-Effects Models with Julia"

This repository uses [Quarto](https://quarto.org) with the Julia code execution supplied by [QuartoNotebookRunner.jl](https://github.com/PumasAI/QuartoNotebookRunner.jl/), which requires Quarto 1.5+.

As of early May 2024, Quarto 1.5 is only available as a preview release, which you'll need to [download from GitHub](https://github.com/quarto-dev/quarto-cli/releases). Under each release's "Assets", you can find platform-specific installers.

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
