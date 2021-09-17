function dyestufftable()
    res = @combine(
        groupby(DataFrame(MixedModels.dataset(:dyestuff)), :batch),
        :mean_yield = mean(:yield),
        :n = length(:yield),
    )
end

function dyestuffdataplot()
    CairoMakie.activate!()
    dyestuff = DataFrame(MixedModels.dataset(:dyestuff))
    myld = @combine(groupby(dyestuff, :batch), :mean_yield=mean(:yield))
    ord = sortperm(myld.mean_yield)
    obatch = CategoricalArray(dyestuff.batch; levels=myld.batch[ord], ordered=true)
    axis = (;
        xlabel="Yield of dyestuff [g]",
        ylabel="Batch of intermediate product",
        yticks=(1:6, levels(obatch)),
    )
    scatter(Array(dyestuff.yield), obatch.refs; color=(:blue, 0.3), axis)
    lines!(myld.mean_yield[ord], 1:6)
    filename = "dyestuff_data"
    caption="Yield of dyestuff by batch.  The line joins the mean yields by batch."
    label = "dyestuffdata"
    Options(current_figure();  filename, caption, label)
end

function dyestuff2dataplot()
    dyestuff2 = select(    # convert columns to types that Makie knows how to deal with
        DataFrame(MixedModels.dataset(:dyestuff2)),
        :yield => Array,
        :batch => PooledArray;
        renamecols=false,
    )
    gdf = groupby(dyestuff2, :batch)
    mnyld = @combine(gdf, :mean_yield=mean(:yield))
    perm = sortperm(mnyld.mean_yield)
    iperm = invperm(perm)
    scatter(
        dyestuff2.yield,
        iperm[dyestuff2.batch.refs] .+ randn(30) .* 0.05,
        axis=(;
            xlabel="Simulated yield (dimensionless)",
            ylabel="Batch in simulation",
            yticks=(1:6, levels(dyestuff2.batch)[perm]),
        ),
    )
    lines!(mnyld.mean_yield[perm], 1:6)
    Options(current_figure();
        caption="Artificial data of yield by batch.  The line joins the mean yields.",
        label="dyestuff2data",
    )
end
