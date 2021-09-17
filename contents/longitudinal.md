# Models for Longitudinal Data {#chap:longitudinal}
============================

Longitudinal data consist of repeated measurements on the same subject (or some other "experimental unit") taken over time.
Generally we wish to characterize the time trends within subjects and between subjects.
The data will always include the response, the time covariate and the indicator of the subject on which the measurement has been made.
If other covariates are recorded, say whether the subject is in the treatment group or the control group, we may wish to relate the within- and between-subject trends to such covariates.

In this chapter we introduce graphical and statistical techniques for the analysis of longitudinal data by applying them to a simple example.

The *sleepstudy* Data {#sec:sleep}
---------------------

@belenky03:_patter report on a study of the effects of sleep deprivation on reaction time for a number of subjects chosen from a population of long-distance truck drivers.
These subjects were divided into groups that were allowed only a limited amount of sleep each night.
We consider here the group of 18 subjects who were restricted to three hours of sleep per night for the first ten days of the trial.
Each subject's reaction time was measured several times on each day of the trial.

```jl
sco("sleepstudy = MixedModels.dataset(:sleepstudy)")
```

```jl
s = """
sleepstudy = DataFrame(sleepstudy)
describe(sleepstudy)
"""
sco(s)
```

In this data frame, the response variable, `reaction`, is the average of the reaction time measurements on a given subject for a given day.
The two covariates are `days`, the number of days of sleep deprivation, and `subj`, the identifier of the subject on which the observation was made.

As recommended for any statistical analysis, we begin by plotting the data.
The most important relationship to plot for longitudinal data on multiple subjects is the trend of the response over time by subject, as shown in @fig:sleepxyplot.
This plot, in which the data for different subjects are shown in separate panels with the axes held constant for all the panels, allows for examination of the time-trends within subjects and for comparison of these patterns between subjects.
Through the use of small panels in a repeating pattern @fig:sleepxyplot conveys a great deal of information, the individual time trends for 18 subjects over 10 days --- a total of 180 points --- without being overly cluttered.

### Characteristics of the *sleepstudy* data {#sec:DataPlotChar}

The principles of "Trellis graphics", developed by Bill Cleveland and his coworkers at Bell Labs and implemented in the `lattice` package for [R](https://R-project.org) by Deepayan Sarkar, have been incorporated in this plot.
As stated above, all the panels have the same vertical and horizontal scales, allowing us to evaluate the pattern over time for each subject and also to compare patterns between subjects.
The line drawn in each panel is a simple least squares line fit to the data in that panel only.
It is provided to enhance our ability to discern patterns in both the slope (the typical change in reaction time per day of sleep deprivation for that particular subject) and the intercept (the average response time for the subject when on their usual sleep pattern).

The aspect ratio of the panels (ratio of the height to the width) has
been chosen, according to an algorithm described in @cleveland93:_visual_data, to facilitate comparison of slopes.
The effect of choosing the aspect ratio in this way is to have the slopes of the lines on the page distributed around $\pm 45^\circ$, thereby making it easier to detect systematic changes in slopes.

The panels have been ordered (from left to right starting at the bottom row) by increasing intercept.
Because the subject identifiers, shown in the strip above each panel, are unrelated to the response it would not be helpful to use the default ordering of the panels, which is by increasing subject number.
If we did so our perception of patterns in the data would be confused by the, essentially random, ordering of the panels.
Instead we use a characteristic of the data to determine the ordering of the panels, thereby enhancing our ability to compare across panels.
For example, a question of interest to the experimenters is whether a subject's rate of change in reaction time is related to the subject's initial reaction time.
If this were the case we would expect that the slopes would show an increasing trend (or, less likely, a decreasing trend) in the left to right, bottom to top ordering.

There is little evidence in @fig:sleepxyplot of such a systematic relationship between the subject's initial reaction time and their rate of change in reaction time per day of sleep deprivation.
We do see that for each subject, except `S335`, reaction time increases, more-or-less linearly, with days of sleep deprivation.
However, there is considerable variation both in the initial reaction time and in the daily rate of increase in reaction time.
We can also see that these data are balanced, both with respect to the number of observations on each subject, and with respect to the times at which these observations were taken.

In cases like this where there are several observations (10) per subject and a relatively simple within-subject pattern (more-or-less linear) we may want to examine coefficients from within-subject fixed-effects fits.
However, because the subjects constitute a sample from the population of interest and we wish to drawn conclusions about typical patterns in the population and the subject-to-subject variability of these patterns, we will eventually want to fit mixed models and we begin by doing so.
In we will compare estimates from a mixed-effects model with those from the
within-subject fixed-effects fits.

## Mixed-effects models For the *sleepstudy* Data {#sec:SleepMixed}

Based on our preliminary graphical exploration of these data, we fit a mixed-effects model with two fixed-effects parameters, the intercept and slope of the linear time trend for the population, and two random effects for each subject.
The random effects for a particular subject are the deviations in intercept and slope of that subject's time trend from the population values.

We will fit two linear mixed models to these data.
One model, `m11`, allows for correlation (in the unconditional distribution) of the random effects for the same subject.
That is, we allow for the possibility that, for example, subjects with higher initial reaction times may, on average, be more strongly affected by sleep deprivation.
The second model. `m12`, provides independent (again, in the unconditional distribution) random effects for intercept and slope for each subject.

### A Model With Correlated Random Effects  {#sec:correlatedre}

The first model is fit as

```jl
s = """
m11 = fit(
    MixedModel,
    @formula(reaction ~ 1 + days + (1+days|subj)),
    MixedModels.dataset(:sleepstudy);
    thin=1,
)
"""
sco(s)
```

From the display we see that this model incorporates both an intercept and a slope (with respect to ) in the fixed effects and in the random effects.
Extracting the conditional modes of the random effects

```jl
sco("first(m11.b)")
```

confirms that these are *vector-valued* random effects.
There are a total of $q=36$ random effects, two for each of the 18 subjects.

The random effects section of the model display,

```jl
sco("VarCorr(m11)")
```

indicates that there will be a random effect for the intercept and a
random effect for the slope with respect to at each level of and,
furthermore, the unconditional distribution of these random effects,
$\mathcal{B}\sim\mathcal{N}(\mathbf{0},\Sigma)$, allows for correlation of the
random effects for the same subject.

We can confirm the potential for correlation of random effects within
subject in the images of $\Lambda$, $\Sigma$ and $\mathbf{L}$ for this model

```jl
sco("only(m11.reterms).λ")
```

The matrix $\Lambda$ has 18 triangular blocks of size 2 along the diagonal, generating 18 square, symmetric blocks of size 2 along the diagonal of $\Sigma$.
The 18 symmetric blocks on the diagonal of $\Sigma$ are identical.
Overall we estimate two standard deviations and a correlation for a vector-valued random effect of size 2, as shown in the model summary.

Often the variances and the covariance of random effects are quoted, rather than the standard deviations and the correlation shown here.
We have already seen that the variance of a random effect is a poor scale on which to quote the estimate because confidence intervals on the variance are so badly skewed.
It is more sensible to assess the estimates of the standard deviations of random effects or, possibly, the logarithms of the standard deviations if we can be confident that 0 is outside the region of interest.
We do display the estimates of the variances of the random effects but mostly so that the user can compare these estimates to those from other software or for cases where an estimate of a variance is expected (sometimes even required) to be given when reporting a mixed model fit.

We do not quote estimates of covariances of vector-valued random effects because the covariance is a difficult scale to interpret, whereas a correlation has a fixed scale.
A correlation must be between $-1$ and $1$, allowing us to conclude that a correlation estimate close to those extremes indicates that $\Sigma$ is close to singular and the model is not well formulated.

The estimates of the fixed effects parameters are $\widehat{\beta}=(251.41,10.467)\trans$.
These represent a typical initial reaction time (i.e. without sleep deprivation) in the population of about 250 milliseconds, or 1/4 sec., and a typical increase in reaction time of a little more than 10 milliseconds per day of sleep deprivation.

The estimated subject-to-subject variation in the intercept corresponds to a standard deviation of about 25 ms.
A 95% prediction interval on this random variable would be approximately $\pm 50$ ms.
Combining this range with a population estimated intercept of 250 ms. indicates that we should not be surprised by intercepts as low as 200 ms. or as high as 300 ms.
This range is consistent with the reference lines shown in @fig:sleepxyplot .

Similarly, the estimated subject-to-subject variation in the slope corresponds to a standard deviation of about 5.7 ms./day so we would not be surprised by slopes as low as $10.5 - 2\cdot 5.7=-0.9$ ms./day or as high as $10.5 + 2\cdot 5.7=21.9$ ms./day.
Again, the conclusions from these rough, "back of the envelope" calculations are consistent with our observations from @fig:sleepxyplot .

The estimated residual standard deviation is about 25 ms. leading us to expect a scatter around the fitted lines for each subject of up to $\pm 50$ ms.
From @fig:sleepxyplot we can see that some subjects (`S309`, `S372`, and `S337`) appear to have less variation than $\pm 50$ ms. about their within-subject fit but others (`S308`, `S332`, and `S331`) may have more.

Finally, we see the estimated within-subject correlation of the random effect for the intercept and the random effect for the slope is very low, $0.081$, confirming our impression that there is little evidence of a systematic relationship between these quantities.
In other words, observing a subject's initial reaction time does not give us much information for predicting whether their reaction time will be strongly affected by each day of sleep deprivation or not.
It seems reasonable that we could get nearly as good a fit from a model that does not allow for correlation, which we describe next.

### A Model With uncorrelated random effects {#sec:uncorrelatedre}

In a model with uncorrelated random effects we have $\mathcal{B}\sim\mathcal{N}(\mathbf{0},\Sigma)$ where $\Sigma$ is diagonal.
We have seen models like this in previous chapters but those models had simple scalar random effects for all the grouping factors.
Here we want to have a simple scalar random effect for and a random effect for the slope with respect to `days`, also indexed by `subj`.
We accomplish this by specifying two random-effects terms.
The first, `(1|subj)`, is a simple scalar term.
The second has `days` on the left hand side of the vertical bar.

It may seem that the model formula we want should be
```
reaction ~ 1 + days + (1|subj) + (days|subj)
```
but it is not.
Because the intercept is implicit in linear models, the second random effects term is equivalent to `(1+days|subj)` and will, by itself, produce
correlated, vector-valued random effects.

We must suppress the implicit intercept in the second random-effects term, which we do by writing it as `(0+days|subj)`, read as "no intercept and `days` by `subj`".
Using the first form we have

```jl
s = """
m12 = fit(
    MixedModel,
    @formula(reaction ~ 1 + days + (1|subj) + (0+days|subj)),
    sleepstudy;
    thin=1,
)
"""
sco(s)
```

As in model `m11`, there are two random effects for each subject

```jl
sco("only(m12.b)")
```

but no correlation has been estimated

```jl
sco("VarCorr(m12)")
```

Images of the matrices $\Lambda$, $\Sigma$ and $\mathbf{L}$ show that $\Sigma$ is indeed diagonal.

Images of $\mathbf{Z}'$ for these two models

```jl
sco("only(m11.reterms).adjA")
```

```jl
sco("only(m12.reterms).adjA")
```

shows that the columns of $\mathbf{Z}$ (rows of
$\mathbf{Z}\trans$) from one model are the same those from the other model.

### Generating $\mathbf{Z}$ and $\Lambda$ From random-effects terms {#sec:GeneratingZLambda}

The non-zero values in the model matrix $\mathbf{Z}$ for model are the same as those for model but the columns are in a different order.
Pairs of columns associated with the same level of the grouping factor are adjacent.
One way to think of the process of generating these columns is to extend the idea of an interaction between a single covariate and the grouping factor to generating an "interaction" of a model matrix and the levels of the grouping factor.
In other words, we begin with the two columns of the model matrix for the expression and the 18 columns of indicators for the factor.
The result will have 36 columns that we regard as 18 adjacent pairs.
The values within each of these pairs of columns are the values of the columns, when the indicator is 1, otherwise they are zero.

We can now describe the general process of creating the model matrix, $\mathbf{Z}$, and the relative covariance factor, $\Lambda$ from the random-effects terms in the model formula.
Each random-effects term is of the form `(A|F)`.
The expression `A` is evaluated as a linear model formula, producing a model matrix with $s$ columns.
The expression `F` is evaluated as a factor.
Let $k$ be the number of levels in this factor, after eliminating unused levels, if any.
The $i$th term generates $s_ik_i$ columns in the model matrix, $\mathbf{Z}$, and a diagonal block of size $s_ik_i$ in the relative covariance factor, $\Lambda$.
The $s_ik_i$ columns in $\mathbf{Z}$ have the pattern of the interaction of the $s_i$ columns from the $i$th with the $k_i$ indicator columns for the $i$th grouping factor `F`.
The diagonal block in $\Lambda$ is itself block diagonal, consisting of $k_i$ blocks, each a lower triangular matrix of size $s_i$.
In fact, these inner blocks are repetitions of the same lower triangular $s_i\times s_i$ matrix.
The $i$ term contributes $s_i(s_i+1)/2$ elements to the variance-component parameter, $\vec\theta$, and these are the elements in the lower triangle of this $s_i\times s_i$ template matrix.

Note that when counting the columns in a model matrix we must take into account the implicit intercept term.
For example, we could write the formula for model as
```
reaction ~ days + (days|subj)
```
realizing that the linear model expression, `days`, actually generates two columns because of the implicit intercept.

Whether or not to include an explicit intercept term in a model formula is a matter of personal taste.
Many people prefer to write the intercept explicitly so as to emphasize the relationship between terms in the formula and coefficients or random effects in the model.
Others omit these implicit terms so as to economize on the amount of typing required.
Either approach can be used.
The important point to remember is that the intercept must be explicitly suppressed when you don't want it in a term.

Also, the intercept term must be explicit when it is the only term in the expression.
That is, a simple, scalar random-effects term must be written as `(1|F)` because a term like `(|F)` is not syntactically correct.
However, we can omit the intercept from the fixed-effects part of the model formula if we have any random-effects terms.
That is, we could write the formula for model in Chap. [\[chap:ExamLMM\]](#chap:ExamLMM){reference-type="ref" reference="chap:ExamLMM"} as

or even

although omitting the parentheses around a random-effects term is risky.
Because of operator precedence, the vertical bar operator, , takes essentially everything in the expression to the left of it as its first operand.
It is advisable always to enclose such terms in parentheses so the scope of the operands to the operator is clearly defined.

### Comparing Models *m11* and *m12* {#sec:comparingfm06fm07}

Returning to models `m11` and `m12` for the data, it is easy to see that these are nested models because `m11` is reduced to `m12` by constraining the within-group correlation of random effects to be zero (which is equivalent to constraining the element below the diagonal in the $2\times 2$ lower triangular $\lambda$ to be zero).

We can use a likelihood ratio test to compare these fitted models.

```jl
sco("MixedModels.likelihoodratiotest(m12, m11)")
```

The value of the $\chi^2$ statistic, $0.0639$, is very small, corresponding to a p-value of $0.80$ and indicating that the extra parameter in model relative to does not produce a significantly better fit.
By the principal of parsimony we prefer the reduced model, `m12`.

This conclusion is consistent with the visual impression provided by @fig:sleepxyplot.
There does not appear to be a strong relationship between a subject's initial reaction time and the extent to which his or her reaction time is affected by sleep deprivation.

In this likelihood ratio test the value of the parameter being tested, a correlation of zero, is not on the boundary of the parameter space.
We can be confident that the p-value from the LRT adequately reflects the underlying situation.

## Assessing the Precision of the Parameter Estimates {#sec:assess-prec-param}

Plots of the profile $\zeta$ for the parameters in model (Fig. [\[fig:pr07plt\]](#fig:pr07plt){reference-type="ref" reference="fig:pr07plt"}) show that confidence intervals on $\sigma_1$ and $\sigma_2$ will be slightly skewed; those for $\log(\sigma)$ will be symmetric and well-approximated by methods based on quantiles of the standard normal distribution and those for the fixed-effects parameters, $\beta_1$ and $\beta_2$ will be symmetric and slightly over-dispersed relative to the standard normal.
For example, the 95% profile-based confidence intervals are

The profile pairs plot (Fig. [\[fig:pr07pairs\]](#fig:pr07pairs){reference-type="ref" reference="fig:pr07pairs"})

shows, for the most part, the usual patterns.
First, consider the panels below the diagonal, which are on the $(\zeta_i,\zeta_j)$ scales.
The $\zeta$ pairs for $\log(\sigma)$ and $\beta_0$, in the $(4,3)$ panel, and for $\log(\sigma)$ and $\beta_1$, in the $(5,3)$ panel, show the ideal pattern.
The profile traces are straight and orthogonal, producing interpolated contours on the $\zeta$ scale that are concentric circles centered at the origin.
When mapped back to the scales of $\log(\sigma)$ and $\beta_0$ or $\beta_1$, in panels $(3,4)$ and $(3,5)$, these circles become slightly distorted, but this is only due to the moderate nonlinearity in the profile $\zeta$ plots for these parameters.

Examining the profile traces on the $\zeta$ scale for $\log(\sigma)$ versus $\sigma_1$, the $(3,1)$ panel, or versus $\sigma_2$, the $(3,2)$ panel, and for $\sigma_1$ versus $\sigma_2$, the $(2,1)$ panel, we see that close to the estimate the traces are orthogonal but as one variance component becomes small there is usually an increase in the others.
In some sense the total variability in the response will be partitioned across the contribution of the fixed effects and the variance components.
In each of these panels the fixed-effects parameters are at their optimal values, conditional on the values of the variance components, and the variance components must compensate for each other.
If one is made smaller then the others must become larger to compensate.

The patterns in the $(4,1)$ panel ($\sigma_1$ versus $\beta_0$, on the $\zeta$ scale) and the $(5,2)$ panel ($\sigma_2$ versus $\beta_1$, on the $\zeta$ scale) are what we have come to expect.
As the fixed-effects parameter is moved from its estimate, the standard deviation of the corresponding random effect increases to compensate.
The $(5,1)$ and $(4,2)$ panels show that changing the value of a fixed effect doesn't change the estimate of the standard deviation of the random effects corresponding to the other fixed effect, which makes sense although the perfect orthogonality shown here will probably not be exhibited in models fit to unbalanced data.

In some ways the most interesting panels are those for the pair of fixed-effects parameters: $(5,4)$ on the $\zeta$ scale and $(4,5)$ on the original scale.
The traces are not orthogonal.
In fact the slopes of the traces at the origin of the $(5,4)$ ($\zeta$ scale) panel are the correlation of the fixed-effects estimators ($-0.194$ for this model) and its inverse.
However, as we move away from the origin on one of the traces in the $(5,4)$ panel it curves back toward the horizontal axis (for the horizontal trace) or the vertical axis (for the vertical trace).
In the $\zeta$ scale the individual contours are still concentric ellipses but their eccentricity changes from contour to contour.
The innermost contours have greater eccentricity than the outermost contours.
That is, the outermost contours are more like circles than are the innermost contours.

In a fixed-effects model the shapes of projections of deviance contours onto pairs of fixed-effects parameters are consistent.
In a fixed-effects model the profile traces in the original scale will always be straight lines.
For mixed models these traces can fail to be linear, as we see here, contradicting the widely-held belief that inferences for the fixed-effects parameters in linear mixed models, based on T or F distributions with suitably adjusted degrees of freedom, will be completely accurate.
The actual patterns of deviance contours are more complex than that.

## Examining the Random Effects and Predictions {#sec:fm07re}

The result of applying to fitted linear mixed model is a list of data frames.
The components of the list correspond to the grouping factors in the random-effects terms, not to the terms themselves.
Model `m11` is the first model we have fit with more than one term for the same grouping factor where we can see the combination of random effects from more than one term.

The method for objects produces one plot for each grouping factor.
For scalar random effects the plot is a normal probability plot.
For two-dimensional random effects, including the case of two scalar terms for the same grouping factor, as in this model, the plot is a scatterplot.
For three or more random effects per level of the grouping factor, the plot is a scatterplot matrix.
The left hand panel in Fig. [\[fig:ranefcoeffm07\]](#fig:ranefcoeffm07){reference-type="ref" reference="fig:ranefcoeffm07"} was created with .

The method for a fitted model combines the fixed-effects estimates and the conditional modes of the random effects, whenever the column names of the random effects correspond to the names of coefficients.
For model the fixed-effects coefficients are and and the columns of the random effects match these names.
Thus we can calculate some kind of per-subject "estimates" of the slope and intercept and plot them, as in the right hand panel of Fig. [\[fig:ranefcoeffm07\]](#fig:ranefcoeffm07){reference-type="ref"
reference="fig:ranefcoeffm07"}.
By comparing the two panels in Fig. [\[fig:ranefcoeffm07\]](#fig:ranefcoeffm07){reference-type="ref" reference="fig:ranefcoeffm07"} we can see that the result of the method is simply the conditional modes of the random effects shifted by the coefficient estimates.

It is not entirely clear how we should interpret these values.
They are a combination of parameter estimates with the modal values of random variables and, as such, are in a type of "no man's land" in the probability model.
(In the Bayesian approach [@box73:_bayes_infer_statis_analy] to inference, however, both the parameters and the random effects are random variables and the interpretation of these values is straightforward.)
Despite the difficulties of interpretation in the probability model, these values are of interest because they determine the fitted response for each subject.

Because responses for each individual are recorded on each of ten days we can determine the within-subject estimates of the slope and intercept (that is, the slope and intercept of each of the lines in @fig:sleepxyplot.
In Fig. [\[fig:shrinkage\]](#fig:shrinkage){reference-type="ref" reference="fig:shrinkage"} we compare the within-subject least squares estimates to the per-subject slope and intercept calculated from model .
We see that, in general, the per-subject slopes and intercepts from the mixed-effects model are closer to the population estimates than are the within-subject least squares estimates.
This pattern is sometimes described as a *shrinkage* of coefficients toward the population values.

The term "shrinkage" may have negative connotations.
John Tukey chose to characterize this process in terms of the estimates for individual subjects "borrowing strength" from each other.
This is a fundamental difference in the models underlying mixed-effects models versus strictly fixed-effects models.
In a mixed-effects model we assume that the levels of a grouping factor are a selection from a population and, as a result, can be expected to share characteristics to some degree.
Consequently, the predictions from a mixed-effects model are attenuated relative to those from strictly fixed-effects models.

The predictions from model and from the within-subject least squares fits for each subject are shown in Fig. [\[fig:shrinkfit\]](#fig:shrinkfit){reference-type="ref" reference="fig:shrinkfit"}.

It may seem that the shrinkage from the per-subject estimates toward the population estimates depends only on how far the per-subject estimates (solid lines) are from the population estimates (dot-dashed lines).
However, careful examination of this figure shows that there is more at work here than a simple shrinkage toward the population estimates proportional to the distance of the per-subject estimates from the population estimates.

It is true that the mixed model estimates for a particular subject are "between" the within-subject estimates and the population estimates, in the sense that the arrows in Fig. [\[fig:shrinkage\]](#fig:shrinkage){reference-type="ref" reference="fig:shrinkage"} all point somewhat in the direction of the population estimate.
However, the extent of the attenuation of the within-subject estimates toward the population estimates is not simply related to the distance between those two sets of estimates.
Consider the two panels, labeled 330 and 337, at the top right of Fig. [\[fig:shrinkfit\]](#fig:shrinkfit){reference-type="ref" reference="fig:shrinkfit"}.
The within-subject estimates for 337 are quite unlike the population estimates but the mixed-model estimates are very close to these within-subject estimates.
That is, the solid line and the dashed line in that panel are nearly coincident and both are a considerable distance from the dot-dashed line.
For subject 330, however, the dashed line is more-or-less an average of the solid line and the dot-dashed line, even though the solid and dot-dashed lines are not nearly as far apart as they are for subject 337.

The difference between these two cases is that the within-subject estimates for 337 are very well determined.
Even though this subject had an unusually large intercept and slope, the overall pattern of the responses is very close to a straight line.
In contrast, the overall pattern for 330 is not close to a straight line so the within-subject coefficients are not well determined.
The multiple $R^2$ for the solid line in the 337 panel is $93.3\%$ but in the 330 panel it is only $15.8\%$.
The mixed model can pull the predictions in the 330 panel, where the data are quite noisy, closer to the population line without increasing the residual sum of squares substantially.
When the within-subject coefficients are precisely estimated, as in the 337 panel, very little shrinkage takes place.

We see from Fig. [\[fig:shrinkfit\]](#fig:shrinkfit){reference-type="ref" reference="fig:shrinkfit"} that the mixed-effects model smooths out the  between-subject differences in the predictions by bringing them closer to a common set of predictions, but not at the expense of dramatically increasing the sum of squared residuals.
That is, the predictions are determined so as to balance fidelity to the data, measured by the residual sum of squares, with simplicity of the model.
The simplest model would use the same prediction in each panel (the dot-dashed line) and the most complex model, based on linear relationships in each panel, would correspond to the solid lines. The dashed lines are between these
two extremes.
We will return to this view of the predictions from mixed models balancing complexity versus fidelity in , where we make the mathematical nature of this balance explicit.

We should also examine the prediction intervals on the random effects (Fig. [\[fig:caterpillar\]](#fig:caterpillar){reference-type="ref" reference="fig:caterpillar"}) where we see that many prediction intervals overlap zero but there are several that do not.
In this plot the subjects are ordered from bottom to top according to increasing conditional mode of the random effect for `(1|subj)`.
The resulting pattern in the conditional modes of the random effect for `(0+days|subj)` reinforces our conclusion that the model `m12`, which does not allow for correlation of the random effects for `(1|subj)` and `(0+days|subj)`, is suitable.

## Chapter Summary
