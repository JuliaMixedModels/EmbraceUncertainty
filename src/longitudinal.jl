function egspaghettiplot(egdf)
    egplt = data(egdf) * mapping(:time => "Age (yr)", :resp => "Ramus bone length (mm)")
	draw(egplt * mapping(color=:Subj) * (visual(Scatter) + visual(Lines)))
    return Options(
        current_figure();
        caption="Length of ramus bone versus age for a sample of 20 boys.",
        label="egspaghetti",
    )
end

function eglayout!(egdf)
    if "orderedsubj" ∉ names(egdf)
        egdf.orderedsubj = CategoricalArray(
    		egdf.Subj; 
	    	levels = @sort(@subset(egdf, :time == 8.0), :resp).Subj,
	    )
    end
    egplt = data(egdf) * mapping(:time => "Age (yr)", :resp => "Ramus bone length (mm)")
    egplt2 = egplt * mapping(layout=:orderedsubj) * visual(ScatterLines)
    draw(egplt2; axis=(height=300, width=120))
    return Options(current_figure();
        caption="Length of ramus bone versus age for a sample of 20 boys.  The panels are ordered rowwise, starting at the top left, by increasing initial bone length",
        label="eglayout",
    )
end

function bxctrlspaghetti(bxdf, bxaxes)
    draw(
        data(@subset(bxdf, :Group == "Control")) * 
        bxaxes * mapping(; color=:Subj) *
        (visual(Scatter) + visual(Lines))
    )
    return Options(
        current_figure(),
        caption="Weight (g) of rats in the control group of bxdf versus time in trial (wk).",
        label="bxctrlspaghetti",
    )
end

function bxctrllayout(bxdf, bxaxes)
    df = @subset(bxdf, :Group == "Control")
    df.orderedsubj = CategoricalArray(
        df.Subj;
        levels=sort(@subset(df, iszero(:time)), :resp).Subj,
        ordered=true,
    )
    draw(
        data(df) * bxaxes * mapping(; layout=:orderedsubj) * visual(ScatterLines);
        axis = (height=250,width=150)
    )
    return Options(
        current_figure(),
        caption="Weight (g) of rats in the control group versus time in trial (wk). Panels are ordered by increasing initial weigth.",
        label="bxctrllayout",
    )
end

function bxthiolayout(bxdf, bxaxes)
    df = @subset(bxdf, :Group == "Thioracil")
    df.orderedsubj = CategoricalArray(
        df.Subj;
        levels=sort(@subset(df, iszero(:time)), :resp).Subj,
        ordered=true,
    )
    draw(
        data(df) * bxaxes * mapping(; layout=:orderedsubj) * visual(ScatterLines);
        axis = (height=250,width=150)
    )
    return Options(
        current_figure(),
        caption="Weight (g) of rats in the Thioracil group versus time in trial (wk). Panels are ordered by increasing initial weigth.",
        label="bxthiolayout",
    )
end

function bxthyrlayout(bxdf, bxaxes)
    df = @subset(bxdf, :Group == "Thyroxin")
    df.orderedsubj = CategoricalArray(
        df.Subj;
        levels=sort(@subset(df, iszero(:time)), :resp).Subj,
        ordered=true,
    )
    draw(
        data(df) * bxaxes * mapping(; layout=:orderedsubj) * visual(ScatterLines);
        axis = (height=250,width=150)
    )
    return Options(
        current_figure(),
        caption="Weight (g) of rats in the Thyroxin group versus time in trial (wk). Panels are ordered by increasing initial weigth.",
        label="bxthyrlayout",
    )
end

function fitegm01(egdf)
    return fit(
        MixedModel,
        @formula(resp ~ 1 + time + (1 + time|Subj)),
        egdf;
        contrasts = Dict(
            :Subj => Grouping(),
        ),
        progress = false,
        thin = 1,
    )
end

function fitegm02(egdf)
    return fit(
        MixedModel,
        @formula(resp ~ 1 + time + (1 + time|Subj)),
        egdf;
        contrasts = Dict(
            :Subj => Grouping(),
            :time => Center(8),
        ),
        progress = false,
        thin = 1,
    )
end

function fitegm03(egdf)
    return fit(
        MixedModel,
        @formula(resp ~ 1 + time + (1 + time|Subj)),
        egdf;
        contrasts = Dict(
            :Subj => Grouping(),
            :time => Center(),
        ),
        progress = false,
        thin = 1,
    )
end

function fitegm04(egdf)
    return fit(
        MixedModel,
        @formula(resp ~ 1 + ctime + ctime^2 + (1 + ctime + ctime^2 |Subj)),
        transform(egdf, :time => (x -> x .- 8.75) => :ctime);
        contrasts = Dict(:Subj => Grouping()),
        progress = false,
        thin = 1,
    )
end

function egm01caterpillar(egm01)
    return Options(
        caterpillar(egm01);
        caption="Conditional modes and 95% prediction intervals on the random effects for slope and intercept in model egm01",
        label="egm01caterpillar",
    )
end

function egm02caterpillar(egm02)
    return Options(
        caterpillar(egm02);
        caption="Conditional modes and 95% prediction intervals on the random effects for slope and intercept in model egm02",
        label="egm02caterpillar",
    )
end

function egm02shrinkage(egm02)
    return Options(
        shrinkageplot!(Figure(resolution=(500,500)), egm02);
        caption="Shrinkage plot of the random effects for model egm02.  Blue dots are the conditional means of the random effects at the parameter estimates.  Red dots are the corresponding unconstrained estimates.",
        label="egm02shrinkage",
    )
end

function egm04shrinkage(egm04)
    return Options(
        shrinkageplot!(Figure(resolution=(500,500)), egm04);
        caption="Shrinkage plot of the random effects for model egm04.  Blue dots are the conditional means of the random effects at the parameter estimates.  Red dots are the corresponding unconstrained estimates.",
        label="egm04shrinkage",
    )
end

function fitbxm01(bxdf)
    form = @formula(resp ~ (1 + time + time^2) * Group + (1+time+time^2|Subj))
    return fit(MixedModel, form, bxdf; contrasts=Dict(:Subj => Grouping()))
end

function fitbxm02(bxdf)
    form = @formula(resp ~ 1 + (time + time^2) & Group + (1+time+time^2|Subj))
    return fit(MixedModel, form, bxdf; contrasts=Dict(:Subj => Grouping()))
end

function fitbxm03(bxdf)
    form = @formula(resp ~ 1 + time + time^2 & Group + (1+time+time^2|Subj))
    return fit(MixedModel, form, bxdf; contrasts=Dict(:Subj => Grouping()))
end

function fitbxm04(bxdf)
    form = @formula(resp ~ 1 + time + time^2 & Group + (1+time|Subj))
    return fit(MixedModel, form, bxdf; contrasts=Dict(:Subj => Grouping()))
end

function fitbxm05(bxdf)
    form = @formula(resp ~ 1 + time + time^2 & Group + (1+time|Subj) + (0+time^2|Subj))
    return fit(MixedModel, form, bxdf; contrasts=Dict(:Subj => Grouping()))
end

function bxm03shrinkage(bxm03)
    return Options(
        shrinkageplot!(Figure(resolution=(500,500)), bxm03);
        caption="Shrinkage plot of the random effects for model bxm03.  Blue dots are the conditional means of the random effects at the parameter estimates.  Red dots are the corresponding unconstrained estimates.",
        label="bxm03shrinkage",
    )
end

function bxm03caterpillar(bxm03)
    return Options(
        caterpillar(bxm03);
        caption="Conditional modes and 95% prediction intervals on the random effects for slope and intercept in model bxm03",
        label="bxm03caterpillar",
    )
end

function bxm03plane(bxm03)
    bpts = Point3f.(eachcol(only(bxm03.b)))
    Upts = Point3f.(eachcol(svd(only(bxm03.λ)).U))
    origin = Point3(zeros(Float32, 3))
    xlabel, ylabel, zlabel = only(bxm03.reterms).cnames
    zlabel = "time²"
    perspectiveness = 0.5
    aspect = :data
    f = Figure(resolution=(1200, 500))
    u, v, w = -Upts[2]    # second principle direction flipped to get positive w
    elevation = asin(w)
    azimuth = atan(v, u)
    ax1 = Axis3(f[1,1]; aspect, xlabel, ylabel, zlabel, perspectiveness)
    ax2 = Axis3(f[1,2]; aspect, xlabel, ylabel, zlabel, perspectiveness, elevation, azimuth)
    scatter!(ax1, bpts; marker='∘', markersize=20)
    scatter!(ax2, bpts; marker='∘', markersize=20)
    for p in Upts
        seg = [origin, p]
        lines!(ax1, seg)
        lines!(ax2, seg)
    end
    return Options(
        f;
        caption="Two views of the conditional means of the random effects from model bxm03. The lines from the origin are the principal axes of the unconditional distribution of the random effects.  The panel on the right is looking in along the negative of the second principle axis (red line in left panel).",
        label="bxm03plane",
    )
end

function bxm03rhodens(bxm03pars)
    draw(
        data(@subset(bxm03pars, :type == "ρ")) *
        mapping(:value => "Bootstrap replicates of correlation estimates"; color=:names) *
        AlgebraOfGraphics.density()
    )
    return Options(
        current_figure();
        caption="Kernel density plots of parametric bootstrap estimates of correlation estimates from model bxm03",
        label="bxm03rhodens",
    )
end

function bxm03rhodensatanh(bxm03pars)
    draw(
        data(@transform(@subset(bxm03pars, :type == "ρ"), :z = atanh(clamp(:value, -0.9999, 0.9999)))) *
        mapping(:z => "Fisher's z transformation of correlation estimates"; color=:names) *
        AlgebraOfGraphics.density()
    )
    return Options(
        current_figure();
        caption="Kernel density plots Fisher's z transformation of parametric bootstrap estimates of correlation estimates from model bxm03",
        label="bxm03rhodensatanh",
    )
end

function bxmodfitted(bxm03beta, bxm04beta, bxm05beta)
    times = Float32.((0:256) / 64)
    times² = abs2.(times)
    z257 = zeros(Float32, 257)
    tmat = hcat(
        ones(Float32, 257*3),
        repeat(times, 3),
        vcat(times², z257, z257),
        vcat(z257, times², z257),
        vcat(z257, z257, times²),
    )
    grp = repeat(["Control", "Thioracil", "Thyroxin"], inner=257)
    draw(
        data(
            append!(
                append!(
                    DataFrame(times = tmat[:, 2], wt = tmat * Float32.(bxm03beta), Group = grp, model="bxm03"),
                    DataFrame(times = tmat[:, 2], wt = tmat * Float32.(bxm04beta), Group = grp, model="bxm04"),
                ),
                DataFrame(times = tmat[:, 2], wt = tmat * Float32.(bxm05beta), Group = grp, model="bxm05"),
            )
        ) * mapping(:times => "Time in trial (wk)", :wt => "Weight (gm)", color=:Group, col=:model) *
        visual(Lines);
        axis=(width=600, height=600),
    )
    return Options(
        current_figure(),
        caption="Typical weight curves from models bxm03, bxm04, and bxm05 for each of the three treatment groups.",
        label="bxmodfitted",
    )
end