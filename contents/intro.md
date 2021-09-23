# A Simple, Linear, Mixed-effects Model {#sec:ExamLMM}

In this book we describe the theory behind a type of statistical model called *mixed-effects* models and the practice of fitting and analyzing such models using the [MixedModels](https://github.com/JuliaStats/MixedModels.jl) package for [Julia](https://julialang.org).
These models are used in many different disciplines.
Because the descriptions of the models can vary markedly between disciplines, we begin by describing what mixed-effects models are and by exploring a very simple example of one type of mixed model, the *linear mixed model*.

This simple example allows us to illustrate the use of the `LinearMixedModel` type in the `MixedModels` package for fitting such models and for analyzing the fitted model.
We also describe methods of assessing the precision of the parameter estimates and of visualizing the conditional distribution of the random effects, given the observed data.

## Mixed-effects Models {#sec:memod}

Mixed-effects models, like many other types of statistical models, describe a relationship between a *response* variable and some of the *covariates* that have been measured or observed along with the response.
In mixed-effects models at least one of the covariates is a *categorical* covariate representing experimental or observational "units" in the data set.
In the example from the chemical industry that is given in this chapter, the observational unit is the batch of an intermediate product used in production of a dye.
In medical and social sciences the observational units are often the human or animal subjects in the study.
In agriculture the experimental units may be the plots of land or the specific plants being studied.

In all of these cases the categorical covariate or covariates are observed at a set of discrete *levels*.
We may use numbers, such as subject identifiers, to designate the particular levels that we observed but these numbers are simply labels.
The important characteristic of a categorical covariate is that, at each observed value of the response, the covariate takes on the value of one of a set of distinct levels.

Parameters associated with the particular levels of a covariate are sometimes called the "effects" of the levels.
If the set of possible levels of the covariate is fixed and reproducible we model the covariate using *fixed-effects* parameters.
If the levels that we observed represent a random sample from the set of all possible levels we incorporate *random effects* in the model.

There are two things to notice about this distinction between fixed-effects parameters and random effects.
First, the names are misleading because the distinction between fixed and random is more a property of the levels of the categorical covariate than a property of the effects associated with them.
Secondly, we distinguish between "fixed-effects parameters", which are indeed parameters in the statistical model, and "random effects", which, strictly speaking, are not parameters.
As we will see shortly, random effects are unobserved random variables.

To make the distinction more concrete, suppose that we wish to model the annual reading test scores for students in a school district and that the covariates recorded with the score include a student identifier and the student's gender.
Both of these are categorical covariates.
The levels of the gender covariate, male and female, are fixed.
If we consider data from another school district or we incorporate scores from earlier tests, we will not change those levels.
On the other hand, the students whose scores we observed would generally be regarded as a sample from the set of all possible students whom we could have observed.
Adding more data, either from more school districts or from results on previous or subsequent tests, will increase the number of distinct levels of the student identifier.

*Mixed-effects models* or, more simply, *mixed models* are statistical models that incorporate both fixed-effects parameters and random effects.
Because of the way that we will define random effects, a model with random effects always includes at least one fixed-effects parameter.
Thus, any model with random effects is a mixed model.

We characterize the statistical model in terms of two random variables: a $q$-dimensional vector of random effects represented by the random variable $\mathcal{B}$ and an $n$-dimensional response vector represented by the random variable $\mathcal{Y}$.
(We use upper-case "script" characters to denote random variables.
The corresponding lower-case upright letter denotes a particular value of the random variable.)
We observe the value, $\mathbf{y}$, of $\mathcal{Y}$.
We do not observe the value, $\mathbf{b}$, of $\mathcal{B}$.

When formulating the model we describe the unconditional distribution of $\mathcal{B}$ and the conditional distribution, $(\mathcal{Y}|\mathcal{B}=\mathbf{b})$.
The descriptions of the distributions involve the form of the distribution and the values of certain parameters. 
We use the observed values of the response and the covariates to estimate these
parameters and to make inferences about them.

That's the big picture.
Now let's make this more concrete by describing a particular, versatile class of mixed models called linear mixed models and by studying a simple example of such a model.
First we describe the data in the example.

## The dyestuff and dyestuff2 data {#sec:DyestuffData}

Models with random effects have been in use for a long time.
The first edition of the classic book, *Statistical Methods in Research and Production*, edited by O.L. Davies, was published in 1947 and contained examples of the use of random effects to characterize batch-to-batch variability in chemical processes.
The data from one of these examples are available as the `dyestuff` data in the `MixedModels` package.
In this section we describe and plot these data and introduce a second example, the data, described in @box73:_bayes_infer_statis_analy.

### The dyestuff data {#sec:dyestuff}

The data are described in @davies72:_statis_method_in_resear_and_produc [, Table 6.3, p. 131], the fourth edition of the book mentioned above, as coming from

> an investigation to find out how much the variation from batch to
> batch in the quality of an intermediate product (H-acid) contributes
> to the variation in the yield of the dyestuff (Naphthalene Black 12B)
> made from it. In the experiment six samples of the intermediate,
> representing different batches of works manufacture, were obtained,
> and five preparations of the dyestuff were made in the laboratory from
> each sample. The equivalent yield of each preparation as grams of
> standard colour was determined by dye-trial.

To access these data within `Julia` we must first attach the package to our session using
```julia
using MixedModels
```
The package must be attached before any of the data sets or functions in the package can be used.
If typing this line results in an error report stating that there is no package by this name then you must first install the package.

In what follows, we will assume that the package has been installed and that it has been attached to the session before any of the code shown has been run.

The `MixedModels` package includes several datasets that are used in examples and in tests of the package.
Individual datasets can be loaded with the `dataset` function.
To allow for other packages to incorporate their own `dataset` function without causing name clashes, this function is not exported and must be called with the fully-qualified name `MixedModels.dataset`.
(If this explanation sounds like gibberish to you, the bottom line is that you must use `MixedModels.dataset(:dyestuff)` and not just `dataset(:dyestuff)`.)

```jl
sco("dyestuff = MixedModels.dataset(:dyestuff)")
```

The output indicates that this dataset is a `Table` read from a file in the [Arrow](https://arrow.apache.org) data format.

A `Table` is a simplified tabular form from the [Tables](https://github.com/JuliaData/Tables.jl) package.
It is often convenient to convert the read-only table form to a `DataFrame`

```jl
sc("dyestuff = DataFrame(dyestuff)")
```

The `describe` function in the [DataFrames](https://github.com/JuliaData/DataFrames.jl) package provides a concise description of the structure of the
data,

```jl
sco("describe(dyestuff)", process=without_caption_label)
```

from which we see that it consists of $30$ observations of the `yield`, the response variable, and of the covariate, `batch`, which is a categorical variable whose levels are character strings.

```jl
sco("typeof(dyestuff.batch)")
```

```jl
sco("levels(dyestuff.batch)")
```

If the labels for the factor levels are arbitrary, as they are here, we will use letters instead of numbers for the labels.
That is, we label the batches as `"A"` through `"F"` rather than `1` through `6`.
When the labels are letters it is clear that the variable is categorical.
When the labels are numbers a categorical covariate can be mistaken for a numeric covariate, with unintended consequences.

It is a good practice to apply `describe` to any data frame the first time you work with it and to check carefully that any categorical variables are indeed represented as factors.

The data in a data frame are viewed as a table with columns corresponding to variables and rows to observations.
The functions `first` and `last` the first or last few rows

```jl
sco("first(dyestuff, 7)", process=without_caption_label)
```

```jl
sco("last(dyestuff, 7)", process=without_caption_label)
```
or we could tabulate the data using `DataFrames.groupby` and the `@combine` macro from the [`DataFrameMacros`](https://github.com/jkrumbiegel/DataFrameMacros.jl) package.


```jl
s = """
@combine(groupby(dyestuff, :batch), :mean_yield = mean(:yield), :n = length(:yield))
"""
sc(s)
```

```jl
EU.dyestufftable()
```

Although @tbl:mean_yield does show us an important property of the data, namely that there are exactly $5$ observations on each batch --- a property that we will describe by saying that the data are *balanced* with respect to `batch` --- we usually learn much more about the structure of such data from plots like

```jl
EU.dyestuffdataplot()
```

than we do from numerical summaries.

In @fig:dyestuffdata we can see that there is considerable variability in yield, even for preparations from the same batch, but there is also noticeable batch-to-batch variability.
For example, four of the five preparations from batch F provided lower yields than did any of the preparations from batches C and E.

This plot, and essentially all the other plots in this book, were created using the [Makie](https://makie.juliaplots.io) package (@DanischKrumbiegel2021).

In `yet-to-be-written-appendix` we review some of the principles of data graphics, such as reordering the levels of the factor by increasing mean response, that enhance the informativeness of the plot.
For example, in this plot the levels of `batch` are sorted by increasing mean yield, to make visual comparisons between batches easier, and the vertical positions are *jittered* to avoid overplotting of points.
(Note that the two lowest yields of samples from batch `A` are identical.)
At this point we will concentrate on the information conveyed by the plot and not on how the plot is created.

In @sec:FittingLMMs we will use mixed models to quantify the variability in yield between batches.
For the time being let us just note that the particular batches used in this experiment are a selection or sample from the set of all batches that we wish to consider.
Furthermore, the extent to which one particular batch tends to increase or decrease the mean yield of the process --- in other words, the "effect" of that particular batch on the yield --- is not as interesting to us as is the extent of the variability between batches.
For the purposes of designing, monitoring and controlling a process we want to predict the yield from future batches, taking into account the batch-to-batch variability and the within-batch variability.
Being able to estimate the extent to which a particular batch in the past increased or decreased the yield is not usually an important goal for us.
We will model the effects of the batches as random effects rather than as fixed-effects parameters.

### The dyestuff2 data {#sec:dyestuff2}

The data are simulated data presented in @box73:_bayes_infer_statis_analy [, Table 5.1.4, p. 247] where the authors state

> These data had to be constructed for although examples of this sort
> undoubtedly occur in practice they seem to be rarely published.

The structure and summary

```jl
sco("dyestuff2 = MixedModels.dataset(:dyestuff2)")
```

```jl
s = """
    dyestuff2 = DataFrame(dyestuff2)
    describe(dyestuff2)
    """
sco(s; process=without_caption_label)
```

are intentionally similar to those of the `dyestuff` data.

A data plot in @fig:dyestuffdata:

```jl
EU.dyestuff2dataplot()
```

shows that the batch-to-batch variability in these data is small compared to the within-batch variability.

In some approaches to mixed models it can be difficult to fit models to such data.
Paradoxically, small "variance components" can be more difficult to estimate than large variance components.

The methods we will present are not compromised when estimating small variance components.

## Fitting Linear Mixed Models {#sec:FittingLMMs}

Before we formally define a linear mixed model, let's go ahead and fit models to these data sets using `MixedModels`.
The simplest way to do this is to use the generic `fit` function with arguments describing the type of model to be fit (i.e. `MixedModel`), a *formula* specifying the model and the *data* on which to evaluate the formula.

We will explain the structure of the formula after we have considered an example.

### A model for the dyestuff data {#sec:dyestuffLMM}

We fit a model to the data allowing for an overall level of the `yield` and for an additive random effect for each level of `batch`.

```jl
sco("m1 = fit(MixedModel, @formula(yield ~ 1 + (1|batch)), dyestuff)")
```

The call to `fit` constructs a `LinearMixedModel` object, evaluates the *maximum likelihood* parameter estimates, assigns the results to the name `m1`, and displays a summary of the fitted model.

#### Details of the printed display

The display of the fitted model has four major sections:

1. a description of the model that was fit
2. some statistics characterizing the model fit
3. a summary of properties of the random effects and
4. a summary of the fixed-effects parameter estimates.

We consider each of these sections in turn.

The description section states that this is a linear mixed model in which the parameters have been estimate by maximum likelihood (ML).
The formula argument is displayed for later reference.

The display of a model fit by maximum likelihood (ML) provides several other model-fit statistics such as Akaike's Information Criterion [@saka:ishi:kita:1986], Schwarz's Bayesian Information Criterion [@schw:1978], the log-likelihood at the parameter estimates, and negative twice the log-likelihood, which is the estimation criterion transformed to the scale of the *deviance*.
For linear mixed models we refer to `-2 loglik` as the value of the *objective* because this is the value that is minimized during the optimization phase of fitting the model.
To evaluate the *deviance* we should subtract the value of this criterion at a *saturated* or baseline model but it is not clear how to define such a baseline model in these cases.

However, it is still possible to perform *likelihood ratio tests* of different models fit to the same data using the difference in the minimized objectives, because it is the same as the difference in the deviances.
(Recall that the objective is negative twice the log-likelihood, hence a ratio of likelihoods corresponds to the difference in objectives.)

The third section is the table of estimates of parameters associated with the random effects.
There are two sources of variability in the model we have fit, a batch-to-batch variability in the level of the response and the residual or per-observation variability --- also called the within-batch variability.
The name "residual" is used in statistical modeling to denote the part of the variability that cannot be explained or modeled with the other terms.
It is the variation in the observed data that is "left over" after we have determined the estimates of the parameters in the other parts of the model.

Some of the variability in the response is associated with the fixed-effects terms.
In this model there is only one such term, labeled as the `(Intercept)`.
The name "intercept", which is better suited to models based on straight lines written in a slope/intercept form, should be understood to represent an overall "typical" or mean level of the response in this case.
(In case you are wondering about the parentheses around the name `(Intercept)`, they are included so that you can't accidentally create a variable with a name that conflicts with this name.)
The line labeled `batch` in the random effects table shows that the random effects added to the term, one for each level of the factor `batch`, are modeled as random variables whose unconditional variance is estimated as 1388.33 g$^2$ in the ML fit.
The corresponding standard deviation is 37.36 g.

Note that the last column in the random effects summary table is the estimate of the variability expressed as a standard deviation rather than as a variance.
These values are provided because it is usually easier to visualize standard deviations, which are on the scale of the response, than it is to visualize the magnitude of a variance.
The values in this column are a simple re-expression (the square root) of the estimated variances.
Do not confuse them with the standard errors of the variance estimators, which are not given here.
In `add-section-reference-here` we explain why we do not provide standard errors of variance estimates.

The line labeled `Residual` in this table gives the estimate of the variance of the residuals (also in g$^2$) and its corresponding standard deviation.
The estimated standard deviation of the residuals is 49.5 g.

The last line in the random effects table states the number of observations to which the model was fit and the number of levels of any "grouping factors" for the random effects.
In this case we have a single random effects term, `(1|batch)`, in the model formula and the grouping factor for that term is `batch`.
There will be a total of six random effects, one for each level of `batch`.

The final part of the printed display gives the estimates and standard errors of any fixed-effects parameters in the model.
The only fixed-effects term in the model formula is the `1`, denoting a constant which, as explained above, is labeled as `(Intercept)`.
The estimate of this parameter is 1527.5 g, which happens to be the mean yield across all the data.

### A model for the dyestuff2 data {#sec:Dyestuff2LMM}

Fitting a similar model to the data produces an estimate $\widehat{\sigma}_1=0$.

```jl
sco("m2 = fit(MixedModel, @formula(yield ~ 1 + (1|batch)), dyestuff2)")
```

An estimate of $0$ for $\sigma_1$ does not mean that there is no variation between the groups.
Indeed @fig:dyestuff2data shows that there is some small amount of variability between the groups.
The estimate, $\widehat{\sigma}_1=0$, simply indicates that the level of "between-group" variability is not sufficient to warrant incorporating random effects in the model.

The important point to take away from this example is that we must allow for the estimates of variance components to be zero.
We describe such a model as being degenerate, in the sense that it corresponds to a linear model in which we have removed the random effects associated with `batch`.
Degenerate models can and do occur in practice.
Even when the final fitted model is not degenerate, we must allow for such models when determining the parameter estimates through numerical optimization.

To reiterate, the model corresponds to the linear model because the random effects are inert, in the sense that they have a variance of zero, and hence can be removed.

Notice that the estimate of $\sigma$ from the linear model (called the in the summary) corresponds to the estimate in the REML fit () but not that from the ML fit ().
The fact that the REML estimates of variance components in mixed models generalize the estimate of the variance used in linear models, in the sense that these estimates coincide in the degenerate case, is part of the motivation for the use of the REML criterion for fitting mixed-effects models.

### Further Assessment of the Fitted Models {#sec:furtherassess}

The parameter estimates in a statistical model represent our "best guess" at the unknown values of the model parameters and, as such, are important results in statistical modeling.
However, they are not the whole story.
Statistical models characterize the variability in the data and we must assess the effect of this variability on the parameter estimates and on the precision of predictions made from the model.

In we introduce a method of assessing variability in parameter estimates using the "profiled deviance" and in we show methods of characterizing the conditional distribution of the random effects given the data.
Before we get to these sections, however, we should state in some detail the probability model for linear mixed-effects and establish some definitions and notation.
In particular, before we can discuss profiling the deviance, we should define the deviance.
We do that in the next section.

## The linear mixed-effects probability model {#sec:Probability}

In explaining some of parameter estimates related to the random effects we have used terms such as "unconditional distribution" from the theory of probability.
Before proceeding further we clarify the linear mixed-effects probability model and define several terms and concepts that will be used throughout the book.
Readers who are more interested in practical results than in the statistical theory should feel free to skip this section.

### Definitions and results {#sec:definitions}

In this section we provide some definitions and formulas without derivation and with minimal explanation, so that we can use these terms in what follows.
In @sec:computational we revisit these definitions providing derivations and more explanation.

As mentioned in @sec:memod, a mixed model incorporates two random variables:
$\mathcal{B}$, the $q$-dimensional vector of random effects, and $\mathcal{Y}$, the $n$-dimensional response vector.
In a linear mixed model the unconditional distribution of $\mathcal{B}$ and the conditional distribution, $(\mathcal{Y}|\mathcal{B}=\mathbf{b})$, are both multivariate Gaussian (or "normal") distributions
$$
  \begin{aligned}
    (\mathcal{Y}|\mathcal{B}=\mathbf{b})&\sim\mathcal{N}(\mathbf{X}\mathbf{\beta}+\mathbf{Z}\mathbf{b},\sigma^2\mathbf{I})\\
    \mathcal{B}&\sim\mathcal{N}(\mathbf{0},\Sigma_\theta) .
  \end{aligned}
$$ {#eq:LMMdist}
The *conditional mean* of $\mathcal{Y}$, given $\mathcal{B}=\mathbf{b}$, is the *linear predictor*, $\mathbf{X}\beta+\mathbf{Z}\mathbf{b}$, which depends on the $p$-dimensional *fixed-effects parameter*, $\mathbf{\beta}$, and on $\mathbf{b}$.
The *model matrices*, $\mathbf{X}$ and $\mathbf{Z}$, of dimension $n\times p$ and $n\times q$, respectively, are determined from the formula for the model and the values of covariates.
Although the matrix $\mathbf{Z}$ can be large (i.e. both $n$ and $q$ can be large), it is sparse (i.e. most of the elements in the matrix are zero).

The *relative covariance factor*, $\Lambda_\theta$, is a $q\times q$ matrix, depending on the *variance-component parameter*, $\mathbf{\theta}$, and generating the symmetric $q\times q$ variance-covariance matrix, $\Sigma_\theta$, according to 
$$
\Sigma_\theta=\sigma^2\Lambda_\theta\Lambda_\theta' .
$$ {#eq:relcovfac}
The *spherical random effects*, $\mathcal{U}\sim\mathcal{N}(\mathbf{0},\sigma^2\mathbf{I}_q)$, determine $\mathcal{B}$ according to
$$
\mathcal{B}=\Lambda_\theta\mathcal{U} .
$$

The *penalized residual sum of squares* (PRSS),
$$
r^2(\mathbf{\theta},\mathbf{\beta},\mathbf{u})=
\|\mathbf{y} -\mathbf{X}\mathbf{\beta} -\mathbf{Z}\Lambda_\theta\mathbf{u}\|^2  + \|\mathbf{u}\|^2,
$$
is the sum of the residual sum of squares, measuring fidelity of the model to the data, and a penalty on the size of $\mathbf{u}$, measuring the complexity of the model.
Minimizing $r^2$ with respect to $\mathbf{u}$,
$$
r^2_{\beta,\theta}=\min_{\mathbf{u}}\left\{\|\mathbf{y} -\mathbf{X}\mathbf{\beta} -\mathbf{Z}\Lambda_\theta\mathbf{u}\|^2+\|\mathbf{u}\|^2\right\}
$$
is a direct (i.e. non-iterative) computation during which we calculate the *sparse Cholesky factor*, $\mathbf{L}_\theta$, which is a lower triangular $q\times q$ matrix satisfying
$$
  \mathbf{L}_\theta\mathbf{L}_\theta'=
  \Lambda_\theta'\mathbf{Z}'\mathbf{Z}\Lambda_\theta+\mathbf{I}_q .
$$ {#eq:sparseCholesky1}
where $\mathbf{I}_q$ is the $q\times q$ *identity matrix*.

The objective (negative twice the log-likelihood) for the parameters, given the data, $\mathbf{y}$, is
$$
  d(\mathbf{\theta},\mathbf{\beta},\sigma|\mathbf{y})
  =n\log(2\pi\sigma^2)+\log(|\mathbf{L}_\theta|^2)+\frac{r^2_{\beta,\theta}}{\sigma^2}.
$$ {#eq:LMMdeviance}
where $|\mathbf{L}_\theta|$ denotes the *determinant* of $\mathbf{L}_\theta$.
Because $\mathbf{L}_\theta$ is triangular, its determinant is the product of its diagonal elements.

Because the conditional mean, $\mathbf{\mu}_{\mathcal{Y}|\mathcal{B}=\mathbf{b}}=\mathbf{X}\mathbf{\beta}+\mathbf{Z}\Lambda_\theta\mathbf{u}$, is a linear function of both $\mathbf{\beta}$ and $\mathbf{u}$, minimization of the PRSS with respect to both $\mathbf{\beta}$ and $\mathbf{u}$ to produce
$$
r^2_\theta =\min_{\mathbf{\beta},\mathbf{u}}\left\{\|\mathbf{y} -\mathbf{X}\mathbf{\beta} -\mathbf{Z}\Lambda_\theta\mathbf{u}\|^2+\|\mathbf{u}\|^2\right\}$$
is also a direct calculation.
The values of $\mathbf{u}$ and $\mathbf{\beta}$ that provide this minimum are called, respectively, the *conditional mode*, $\tilde{\mathbf{u}}_\theta$, of the spherical random effects and the conditional estimate, $\widehat{\mathbf{\beta}}_\theta$, of the fixed effects.
At the conditional estimate of the fixed effects the deviance is
$$
  d(\mathbf{\theta},\widehat{\beta}_\theta,\sigma|\mathbf{y})
  =n\log(2\pi\sigma^2)+\log(|\mathbf{L}_\theta|^2)+\frac{r^2_\theta}{\sigma^2}
$$ {#eq:LMMprdev}

Minimizing this expression with respect to $\sigma^2$ produces the conditional estimate
$$
\widehat{\sigma^2}_\theta=\frac{r^2_\theta}{n}
$$
which provides the *profiled objective*,
$$
  \tilde{d}(\mathbf{\theta}|\mathbf{y})=d(\mathbf{\theta},\widehat{\beta}_\theta,\widehat{\sigma}_\theta|\mathbf{y})
  =\log(|\mathbf{L}_\theta|^2)+n\left[1 +
    \log\left(\frac{2\pi r^2_\theta}{n}\right)\right],
$$ {#eq:LMMprofdev}
a function of $\mathbf{\theta}$ alone.

The MLE of $\mathbf{\theta}$, written $\widehat{\mathbf{\theta}}$, is the value that minimizes the profiled objective (#eq:LMMprofdev).
We determine this value by numerical optimization.
In the process of evaluating $\tilde{d}(\widehat{\mathbf{\theta}}|\mathbf{y})$ we determine $\widehat{\mathbf{\beta}}=\widehat{\mathbf{\beta}}_{\widehat\theta}$, $\tilde{\mathbf{u}}_{\widehat{\theta}}$ and $r^2_{\widehat{\theta}}$, from
which we can evaluate $\widehat{\sigma}=\sqrt{r^2_{\widehat{\theta}}/n}$.

The elements of the conditional mode of $\mathcal{B}$, evaluated at the parameter estimates, 
$$
  \tilde{b}_{\widehat{\theta}}=\Lambda_{\widehat{\theta}}\tilde{u}_{\widehat{\theta}}
$$
are sometimes called the *best linear unbiased predictors* or BLUPs of the random effects.
Although it has an appealing acronym, we don't find the term particularly instructive (what is a "linear unbiased predictor" and in what sense are these the "best"?) and prefer the term "conditional mode", which is explained in .

### Matrices and vectors in the fitted model object {#sec:FittedModel}

The optional argument, `thin=1`, in a call to `fit` causes all the values of $\mathbf{\theta}$ and the objective during the progress of the iterative optimization of $\tilde{d}(\mathbf{\theta}|\mathbf{y})$ to be stored in the `optsum` member of the fit.
```jl
sco("""
m1trace = fit(MixedModel, @formula(yield ~ 1 + (1|batch)), dyestuff; thin=1)
m1trace.optsum.fitlog
""")
```
The algorithm converges after 17 function evaluations to a profiled deviance of 327.3270598811301 at $\theta=0.7525807289241839$.
In this model the scalar parameter $\theta$ is the ratio $\sigma_1/\sigma$.

The actual values of many of the matrices and vectors defined above are available as properties of the fitted model object.

In this case the $\Lambda_\theta$ matrix will be a $6\times 6$ diagonal matrix with the diagonal elements all equal to $\theta=0.7525807289241839$.

The Cholesky factor, $\mathbf{L}$, is
```jl
sco("sparseL(m1; full=true)")
```
which consists of three blocks
```jl
sco("BlockDescription(m1)")
```
As we get to larger models we will see that large sparse matrices are displayed as patterns rather than as numerical values.

In this simple model $\Lambda=\theta\mathbf{I}_6$ is a multiple of the identity matrix and the $30\times 6$ model matrix $\mathbf{Z}$, whose transpose is

```jl
sco("Int.(Array(first(m1.reterms)))'")
```

(The conversion to integer, or `Int`, values is to save space when printing.)

Thus $\mathbf{Z}$ is the indicator columns for `batch`.
Because the data are balanced with respect to `batch`, the upper-left block of the Cholesky factor, $\mathbf{L}$, is also a multiple of the identity.

The vector $\mathbf{u}$ is available as a row vector

```jl
sco("first(m1.u)")
```

The vector $\mathbf{\beta}$ and the model matrix $\mathbf{X}$ are available as

```jl
sco("Int.(m1.X')")
```

and

```jl
sco("m1.β")
```

## Assessing the variability of the parameter estimates {#sec:variability}

In this section we show how to create a *profile deviance* object from a fitted linear mixed model and how to use this object to evaluate confidence intervals on the parameters.
We also discuss the construction and interpretation of *profile zeta* plots for the parameters.
In `chapter-that-may-or-may-not-get-written` we discuss the use of the deviance profiles to produce likelihood contours for pairs of parameters.

### Confidence intervals on the parameters {#sec:profdevconf}

The mixed-effects model fit as or has three parameters for which we obtained estimates.
These parameters are $\sigma_1$, the standard deviation of the random effects, $\sigma$, the standard deviation of the residual or "per-observation" noise term and $\beta_0$, the fixed-effects parameter that is labeled as `(Intercept)`.

The function systematically varies the parameters in a model, assessing the best possible fit that can be obtained with one parameter fixed at a specific value and comparing this fit to the *globally optimal fit*, which is the original model fit that allowed all the parameters to vary.
The models are compared according to the change in the deviance, which is the *likelihood ratio test* (LRT) statistic.
We apply a *signed square root* transformation to this statistic and plot the resulting function, called $\zeta$, versus the parameter value.
A $\zeta$ value can be compared to the quantiles of the *standard normal distribution*, $\mathcal{Z}\sim\mathcal{N}(0,1)$.
For example, a 95% profile deviance confidence interval on the parameter consists of those values for which $-1.960 < \zeta < 1.960$.

Because the process of profiling a fitted model, which involves re-fitting the model many times, can be computationally intensive, one should exercise caution with complex models fit to very large data sets.
Because the statistic of interest is a likelihood ratio, the model is re-fit according to the maximum likelihood criterion, even if the original fit is a REML fit.
Thus, there is a slight advantage in starting with an ML fit.

Plots of $\zeta$ versus the parameter being profiled  are obtained with

We will refer to such plots as *profile zeta* plots.
I usually adjust the aspect ratio of the panels in profile zeta plots to, say, and frequently set the layout so the panels form a single row (, in this
case).

The vertical lines in the panels delimit the 50%, 80%, 90%, 95% and 99%
confidence intervals, when these intervals can be calculated.
Numerical values of the endpoints are returned by the extractor.

By default the 95% confidence interval is returned.
The optional argument, , is used to obtain other confidence levels.

Notice that the lower bound on the 99% confidence interval for $\sigma_1$ is not defined.
Also notice that we profile $\log(\sigma)$ instead of $\sigma$, the residual standard deviation.

A plot of $|\zeta|$, the absolute value of $\zeta$, versus the parameter
, obtained by adding the optional argument to the call to , can be more effective for visualizing the confidence intervals.

### Interpreting the profile zeta plot {#sec:interpprofzeta}

A profile zeta plot, such as , shows us the sensitivity of the model fit to changes in the value of particular parameters.
Although this is not quite the same as describing the distribution of an estimator, it is a similar idea and we will use some of the terminology from distributions when describing these plots.
Essentially we view the patterns in the plots as we would those in a normal probability plot of data values or of residuals from a model.

Ideally the profile zeta plot will be close to a straight line over the region of interest, in which case we can perform reliable statistical inference based on the parameter's estimate, its standard error and quantiles of the standard normal distribution.
We describe such a situation as providing a good normal approximation for inference.
The common practice of quoting a parameter estimate and its standard error assumes that this is always the case.

In  the profile zeta plot for $\log(\sigma)$ is reasonably straight so $\log(\sigma)$ has a good normal approximation.
But this does not mean that there is a good normal approximation for $\sigma^2$ or even for $\sigma$.
As shown in  the profile zeta plot for $\log(\sigma)$ is slightly skewed, that for $\sigma$ is moderately skewed and the profile zeta plot for $\sigma^2$ is highly skewed.
Deviance-based confidence intervals on $\sigma^2$ are quite asymmetric, of the form "estimate minus a little, plus a lot".

This should not come as a surprise to anyone who learned in an introductory statistics course that, given a random sample of data assumed to come from a Gaussian distribution, we use a $\chi^2$ distribution, which can be quite skewed, to form a confidence interval on $\sigma^2$.
Yet somehow there is a widespread belief that the distribution of variance estimators in much more complex situations should be well approximated by a normal distribution.
It is nonsensical to believe that.
In most cases summarizing the precision of a variance component estimate by giving an approximate standard error is woefully inadequate.

The pattern in the profile plot for $\beta_0$ is sigmoidal (i.e. an elongated "S"-shape).
The pattern is symmetric about the estimate but curved in such a way that the profile-based confidence intervals are wider than those based on a normal approximation.
We characterize this pattern as symmetric but over-dispersed (relative to a normal distribution).
Again, this pattern is not unexpected. Estimators of the coefficients in a linear model without random effects have a distribution which is a scaled Student's T distribution.
That is, they follow a symmetric distribution that is over-dispersed relative to the normal.

The pattern in the profile zeta plot for $\sigma_1$ is more complex.
 shows the profile zeta plot on the scale of $\log(\sigma_1)$, $\sigma_1$ and $\sigma_1^2$.
Notice that the profile zeta plot for $\log(\sigma_1)$ is very close to linear to the right of the estimate but flattens out on the left.
That is, $\sigma_1$ behaves like $\sigma$ in that its profile zeta plot is more-or-less a straight line on the logarithmic scale, except when $\sigma_1$ is close to zero.
The model loses sensitivity to values of $\sigma_1$ that are close to zero.
If, as in this case, zero is within the "region of interest" then we should expect that the profile zeta plot will flatten out on the left hand side.

Notice that the profile zeta plot of $\sigma_1^2$ in
 is dramatically skewed.
If reporting the estimate, $\widehat{\sigma^2}_1$, and its standard error, as many statistical software packages do, were to be an adequate description of the variability in this estimate then this profile zeta plot should be a straight line.
It's nowhere close to being a straight line in this and in many other model fits, which is why we don't report standard errors for variance estimates.

### Deriving densities from the profile {#sec:profDens}

In the profile zeta plots we show $\zeta$ as a function of a parameter.
We can use the function shown there, which we will call the *profile zeta function*, to generate a corresponding distribution by setting the cumulative distribution function (c.d.f) to be $\Phi(\zeta)$ where $\Phi$ is the c.d.f. of the standard normal distribution.
From this we can derive a density.

This is not quite the same as evaluating the distribution of the estimator of the parameter, which for mixed-effects models can be very difficult, but it gives us a good indication of what the distribution of the estimator would be.

shows the densities corresponding to the profiles in .
We see that the density for $\sigma_1$ is quite skewed.

If we had plotted the densities corresponding to the profiles of the variance components instead, we would get which, of course, just accentuates the skewness in the distribution of these variance components.

## Assessing the random effects {#sec:assessRE}

In @sec:definitions we mentioned that what are sometimes called the BLUPs (or best linear unbiased predictors) of the random effects, $\mathcal{B}$, are the conditional modes evaluated at the parameter estimates, calculated as $\tilde{b}_{\widehat{\theta}}=\Lambda_{\widehat{\theta}}\tilde{u}_{\widehat{\theta}}$.

These values are often considered as some sort of "estimates" of the random effects.
It can be helpful to think of them this way but it can also be misleading.
As we have stated, the random effects are not, strictly speaking, parameters --- they are unobserved random variables.
We don't estimate the random effects in the same sense that we estimate parameters.
Instead, we consider the conditional distribution of $\mathcal{B}$ given the observed data, $(\mathcal{B}|\mathcal{Y}=\mathbf{y})$.

Because the unconditional distribution, $\mathcal{B}\sim\mathcal{N}(\mathbf{0},\Sigma_\theta)$ is continuous, the conditional distribution, $(\mathcal{B}|\mathcal{Y}=\mathbf{y})$ will also be continuous.
In general, the mode of a probability density is the point of maximum density, so the phrase "conditional mode" refers to the point at which this conditional density is maximized.
Because this definition relates to the probability model, the values of the parameters are assumed to be known.
In practice, of course, we don't know the values of the parameters (if we did there would be no purpose in forming the parameter estimates), so we use the estimated values of the parameters to evaluate the conditional modes.

Those who are familiar with the multivariate Gaussian distribution may recognize that, because both $\mathcal{B}$ and $(\mathcal{Y}|\mathcal{B}=\mathbf{b})$ are multivariate Gaussian, $(\mathcal{B}|\mathcal{Y}=\mathbf{y})$ will also be multivariate Gaussian and the conditional mode will also be the conditional mean of $\mathcal{B}$, given $\mathcal{Y}=\mathbf{y}$.
This is the case for a linear mixed model but it does not carry over to other forms of mixed models.
In the general case all we can say about $\tilde{\mathbf{u}}$ or $\tilde{\mathbf{b}}$ is that they maximize a conditional density, which is why we use the term "conditional mode" to describe these values.
We will only use the term "conditional mean" and the symbol, $\mathbf{\mu}$, in reference to $\mathrm{E}(\mathcal{Y}|\mathcal{B}=\mathbf{b})$, which is the conditional mean of $\mathcal{Y}$ given $\mathcal{B}$, and an important part of the formulation of all types of mixed-effects models.

The conditional modes are available as a vector of matrices

```jl
sco("only(m1.b)")
```

In this case the vector consists of a single matrix because there is only one random-effects term, `(1|batch)`, in the model and, hence, only one grouping factor, `batch`, for the random effects.
There is only one row in the matrix because the random-effects term, `(1|batch)`, is a simple, scalar term.

To make this more explicit, random-effects terms in the model formula
are those that contain the vertical bar (`|`) character.
The variable or expression on the right of the `|` is the grouping factor for the random effects generated by this term.
If the expression on the left of the vertical bar is `1`, as it is here, we describe the term as a *simple, scalar, random-effects term*.
The designation "scalar" means there will be exactly one random effect generated for each level of the grouping factor.
A simple, scalar term generates a block of indicator columns --- the indicators for the grouping factor --- in $\mathbf{Z}$.
Because there is only one random-effects term in this model and because that term is a simple, scalar term, the model matrix $\mathbf{Z}$ for this model is the indicator matrix for the levels of `batch`.

In the next chapter we fit models with multiple simple, scalar terms and, in subsequent chapters, we extend random-effects terms beyond simple, scalar terms.
When we have only simple, scalar terms in the model, each term has a unique grouping factor and the elements of the list returned by can be considered as associated with terms or with grouping factors.
In more complex models a particular grouping factor may occur in more than one term.
In such cases the terms associated with the same grouping factor are internally amalgamated into a single term.
Thus internally the random effects are associated with grouping factors, not the terms in the model formula.

Given the data, $\mathbf{y}$, and the parameter estimates, we can evaluate a measure of the dispersion of $(\mathcal{B}|\mathcal{Y}=\mathbf{y})$.
In the case of a linear mixed model, this is the conditional standard deviation, from which we can obtain a prediction interval. 
A plot of these prediction intervals is sometimes called a *caterpillar plot* because it can look like a fuzzy caterpillar when there are many levels of the grouping factor.

```jl
s = """
    CairoMakie.activate!()   # hide
    caterpillar(m1)
    caption = "Caterpillar plot of prediction intervals for m1 random effects" # hide
    label = "caterpillar_m1" # hide
    filename = "caterpillar_m1" # hide
    Options(current_figure(); filename, caption, label) # hide
"""
sco(s)
```

The `caterpillar` function returns a plot with linear spacing of the
levels on the y axis.
An alternative

```jl
s = """
    CairoMakie.activate!()   # hide
    qqcaterpillar(m1)
    caption = "Quantile caterpillar plot of prediction intervals for m1 random effects" # hide
    label = "qqcaterpillar_m1" # hide
    filename = "qqcaterpillar_m1" # hide
    Options(current_figure(); filename, caption, label) # hide
"""
sco(s)
```

returns a plot where the intervals are plotted with vertical spacing corresponding to the quantiles of the standard normal distribution.

The `caterpillar` plot is preferred when there are only a few levels of the grouping factor, as in this case.
When there are hundreds or thousands of random effects the `qqcaterpillar` form is preferred because it focuses attention on the "important few" at the extremes and de-emphasizes the "trivial many" that are close to zero.

## Chapter summary {#sec:ChIntroSummary}

A considerable amount of material has been presented in this chapter, especially considering the word "simple" in its title (it's the model that is simple, not the material).
A summary may be in order.

A mixed-effects model incorporates fixed-effects parameters and random effects, which are unobserved random variables, $\mathcal{B}$.
In a linear mixed model, both the unconditional distribution of $\mathcal{B}$ and the conditional distribution, $(\mathcal{Y}|\mathcal{B}=\mathbf{b})$, are multivariate Gaussian distributions.
Furthermore, this conditional distribution is a spherical Gaussian with mean, $\mathbf{\mu}$, determined by the linear predictor, $\mathbf{Z}\mathbf{b}+\mathbf{X}\mathbf{\beta}$.
That is,
$$
(\mathcal{Y}|\mathcal{B}=\mathbf{b})\sim
  \mathcal{N}(\mathbf{Z}\mathbf{b}+\mathbf{X}\mathbf{\beta}, \sigma^2\mathbf{I}_n) .
$$
The unconditional distribution of $\mathcal{B}$ has mean $\mathbf{0}$ and a parameterized $q\times q$ variance-covariance matrix, $\Sigma_\theta$.

In the models we considered in this chapter, $\Sigma_\theta$, is a simple multiple of the identity matrix, $\mathbf{I}_6$.
This matrix is always a multiple of the identity in models with just one random-effects term that is a simple, scalar term.
The reason for introducing all the machinery that we did is to allow for more general model specifications.

The maximum likelihood estimates of the parameters are obtained by minimizing the deviance.
For linear mixed models we can minimize the profiled deviance, which is a function of $\mathbf{\theta}$ only, thereby considerably simplifying the optimization problem.

To assess the precision of the parameter estimates, we profile the deviance function with respect to each parameter and apply a signed square root transformation to the likelihood ratio test statistic, producing a profile zeta function for each parameter.
These functions provide likelihood-based confidence intervals for the parameters.
Profile zeta plots allow us to visually assess the precision of individual parameters.
Density plots derived from the profile zeta function provide another way of examining the distribution of the estimators of the parameters.

Prediction intervals from the conditional distribution of the random effects, given the observed data, allow us to assess the precision of the random effects.
