function penicillindataplot()
    df, _ = _meanrespfrm(DataFrame(MixedModels.dataset(:penicillin)), :diameter, :plate)
    sort!(df, [:plate, :sample])

    mp = mapping(
        :diameter => "Diameter of inhibition zone [mm]", :plate => "Plate"; color=:sample
    )
    draw(data(df) * mp * visual(ScatterLines; marker='○', markersize=12))
    return Options(
        current_figure();
        caption="Diameter of inhibition zone by plate and sample. Plates are ordered by increasing mean response.",
        label="penicillindot",
    )
end

function fitpnm01()
    return fit(
        MixedModel,
        @formula(diameter ~ 1 + (1 | plate) + (1 | sample)),
        MixedModels.dataset(:penicillin);
        contrasts=Dict(:plate => Grouping(), :sample => Grouping()),
        progress=false,
        thin=1,
    )
end

function pnm01platepred(pnm01)
    caterpillar!(Figure(; resolution=(1000, 600)), ranefinfo(pnm01, :plate))
    return Options(
        current_figure();
        caption="Conditional modes and 95% prediction intervals of random effects for plate in model pnm01",
        filename="pnm01platepred",
        label="pnm01platepred",
    )
end

function pnm01samplepred(pnm01)
    caterpillar!(Figure(; resolution=(1000, 200)), ranefinfo(pnm01, :sample))
    return Options(
        current_figure();
        caption="Conditional modes and 95% prediction intervals of random effects for sample in model pnm01",
        filename="pnm01samplepred",
        label="pnm01samplepred",
    )
end

bssampdens(allpars, type::AbstractString="β") = bssampdens(DataFrame(allpars), type)

function bssampdens(df::DataFrame, type::AbstractString="β")
    selection = Dict("β" => (:names => "Names"), "σ" => (:group => "Group"))
    return draw(
        data(subset(df, :type => x -> x .== type)) *
        mapping(:value => "Bootstrap samples of " * type; color=selection[type]) *
        AlgebraOfGraphics.density(),
    )
end

"""
    objgrid(m::LinearMixedModel, xv::AbstractVector, yv::AbstractVector)

Return `xv`, `yv` and a matrix of objective values for model `m` evaluated on a grid of `θ` values formed by `xv` and `yv`
"""
function objgrid(
    m::LinearMixedModel{T}, xv::AbstractVector{T}, yv::AbstractVector{T}
) where {T}
    θ₀ = m.θ
    θ = copy(θ₀)
    obj = Matrix{T}(undef, (length(xv), length(yv)))
    for (j, y) in enumerate(yv)
        θ[2] = y
        for (i, x) in enumerate(xv)
            θ[1] = x
            obj[i, j] = objective(updateL!(setθ!(m, θ)))
        end
    end
    updateL!(setθ!(m, θ₀))
    return xv, yv, obj
end

function fitpsm01(dat)
    return fit(
        MixedModel,
        @formula(strength ~ 1 + (1 | sample) + (1 | batch)),
        dat;
        contrasts=Dict(:sample => Grouping(), :batch => Grouping()),
        progress=false,
        thin=1,
    )
end

function fitpsm02(dat)
    return fit(
        MixedModel,
        @formula(strength ~ 1 + (1 | sample)),
        dat;
        contrasts=Dict(:sample => Grouping()),
        progress=false,
        thin=1,
    )
end

function iednobshist(insteval)
    return Options(
        draw(
            data(combine(groupby(insteval, :d), nrow => :n)) *
            mapping(:n => "Number of observations") *
            AlgebraOfGraphics.histogram(; bins=410),
        ).figure;
        caption="Histogram of the number of observations per instructor in the insteval data",
        label="iednobshist",
    )
end
