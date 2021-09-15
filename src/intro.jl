function dyestufftable()
    res = @combine(
        groupby(DataFrame(MixedModels.dataset(:dyestuff)), :batch),
        :mean_yield = mean(:yield),
        :n = length(:yield),
    )
    caption = "Mean yield by batch of dyestuff"
    label = "mean_yield"
    Options(res; caption, label)
end

function dyestuffdataplot()
    dyestuff = select(    # convert columns to types that Makie knows how to deal with
        DataFrame(MixedModels.dataset(:dyestuff)),
        :yield => Array,
        :batch => PooledArray;
        renamecols=false,
    )
    gdf = groupby(dyestuff, :batch)
    mnyld = @combine(gdf, :mean_yield=mean(:yield))
    perm = sortperm(mnyld.mean_yield)
    iperm = invperm(perm)
    scatter(
        dyestuff.yield,
        iperm[dyestuff.batch.refs] .+ randn(30) .* 0.05,
        axis=(;
            xlabel="Yield of dyestuff [g]",
            ylabel="Batch of intermediate product",
            yticks=(1:6, levels(dyestuff.batch)[perm]),
        ),
    )
    lines!(mnyld.mean_yield[perm], 1:6)
    Options(current_figure();
        caption="Yield of dyestuff by batch.  The line joins the mean yields.",
        label="dyestuffdata",
    )
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
