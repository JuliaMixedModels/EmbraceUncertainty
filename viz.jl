const hide_progress = true

using AlgebraOfGraphics
using CairoMakie
using DataFrames
if contains(first(Sys.cpu_info()).model, "Intel")
  using MKL
end
using MixedModels
using MixedModelsMakie
using ProgressMeter
using Random

kb07 = DataFrame(MixedModels.dataset("kb07"))
contrasts = Dict(:item => Grouping(),
                 :subj => Grouping(),
                 :spkr => EffectsCoding(),
                 :prec => EffectsCoding(),
                 :load => EffectsCoding())

fm1 = fit(MixedModel, @formula(1000 / rt_raw ~ spkr * prec * load + (1|item) + (1|subj)), kb07; contrasts, progress=true, thin=1)
plot_blups(args...; kwargs...) = plot_blups!(Figure(; resolution=(1000, 600)), args...;kwargs...)
function plot_blups!(f::Makie.FigureLike, m::MixedModel)

    re = ranefinfo(m)
    for (idx, grp) in enumerate(propertynames(re))
      gl = f[1,idx] = GridLayout()
      qqcaterpillar!(gl, re[grp])
      Label(gl[end+1, :], string(grp); font=:bold, tellwidth=false)
    end
    Label(f[0, :], "Conditional Modes", tellwidth=false)
    for idx in 2:length(re)
      ratio = length(re[idx].cnames) / length(re[1].cnames)
      colsize!(f.layout, idx, Auto(ratio))
    end
    return f
end

function eff_dim(m::MixedModel, threshold::Real=0.95)
  re = MixedModels.rePCA(m)
  return NamedTuple{propertynames(re)}(findfirst(>=(threshold), val) for val in values(re))
end

fm_max = fit(MixedModel, @formula(1000 / rt_raw ~ spkr * prec * load + (1 + spkr * prec * load |item) + (1 + spkr * prec * load|subj)), kb07; contrasts, progress=true, thin=1)
fm2 = fit(MixedModel, @formula(1000 / rt_raw ~ spkr * prec * load + (1 + prec + load |item) + (1 + prec + load|subj)), kb07; contrasts, progress=true, thin=1)
fm3 = fit(MixedModel, @formula(1000 / rt_raw ~ spkr * prec * load + (1 + prec + load |item) + (1 + prec + load|subj)), kb07; contrasts, progress=true, thin=1)
fm4 = fit(MixedModel, @formula(1000 / rt_raw ~ spkr * prec * load + (1 + prec + load |item) + (1 + load|subj)), kb07; contrasts, progress=true, thin=1)

function plot_objective(m::MixedModel)
    fitlog = optsum.fitlog
end
