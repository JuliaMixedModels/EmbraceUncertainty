function dyestufftable()
    res = @combine(
        groupby(DataFrame(MixedModels.dataset(:dyestuff)), :batch),
        :mean_yield = mean(:yield),
        :n = length(:yield),
    )
end

function _meanrespfrm(df, resp::Symbol, grps::Symbol; sumryf::Function=mean)
                    # ensure the relevant columns are types that Makie can deal with 
    df = select(df, resp => Array, grps => CategoricalArray; renamecols=false)
                    # create a summary table by mean resp
    sumry = sort!(combine(groupby(df, grps), resp => sumryf => resp), resp)
    glevs = string.(sumry[!, grps])   # group levels in ascending order of mean resp
    levels!(df[!, grps], glevs)
    levels!(sumry[!, grps], glevs)
    return df, sumry
end

function dyestuffdataplot()
    df, sumry = _meanrespfrm(DataFrame(MixedModels.dataset(:dyestuff)), :yield, :batch)

    mp = mapping(:yield => "Yield of dyestuff [g]", :batch => "Batch of intermediate product")
    draw(
        (data(df) * mp * visual(Scatter, marker='○', markersize=12)) +
        (data(sumry) * mp * visual(Lines))
    )
    Options(current_figure();
        caption="Yield of dyestuff by batch.  The line joins the mean yields.",
        label="dyestuffdata",
    )
end

function dyestuff2dataplot()
    df, sumry = _meanrespfrm(DataFrame(MixedModels.dataset(:dyestuff2)), :yield, :batch)
    mp = mapping(:yield => "Simulated yield", :batch => "Batch")
    draw(
        (data(df) * mp * visual(Scatter, marker='○', markersize=12)) +
        (data(sumry) * mp * visual(Lines))
    )
    Options(current_figure();
        caption="Artificial data of yield by batch.  The line joins the mean yields.",
        label="dyestuff2data",
    )
end
