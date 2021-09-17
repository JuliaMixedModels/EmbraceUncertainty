# Generalized Linear Mixed Models for Binary Responses {#sec:GLMMbinomial}

In this chapter we consider mixed-effects models for data sets in which the response is *binary*, representing yes/no or true/false or correct/incorrect responses.

Because the response can only take on one of two values we adapt our models to predict the probability of the positive response.
We retain the concept of the linear predictor, $\mathbf{X}\mathbf{\beta}+\mathbf{Z}\mathbf{b}$, depending on the fixed-effects parameters, $\mathbf{\beta}$, and the random effects, $\mathbf{b}$, with the corresponding model matrices, $\mathbf{X}$ and $\mathbf{Z}$, determining the conditional mean, $\mathbf{\mu}$, of the response, given the random effects, but the relationship is more general than that for the linear mixed model.
The linear predictor, which we shall write as $\mathbf{\eta}$, determines $\mathbf{\mu}$ according to a *link function*, $g$.
For historical reasons it is the function taking an element of $\mathbf{\mu}$ to the corresponding element of $\mathbf{\eta}$ that is called the link.
The transformation in the opposite direction, from $\mathbf{\eta}$ to $\mathbf{\mu}$, is called the *inverse link*.

As described in earlier chapters, models based on a Gaussian distribution for the response vector with its mean determined by the linear predictor are called linear models.
Models incorporating coefficients in a linear predictor but allowing for more general forms of the distribution of the response are called *generalized linear
models*.
When the linear predictor incorporates random effects in addition to the fixed-effects parameters we call them *generalized linear mixed models* (GLMMs) and fit such models with the function.
As in previous chapters, we will begin with an example to help illustrate these ideas.

## Artificial contraception use in regions of Bangladesh {#sec:contraception}

One of the test data sets from the Center for Multilevel Modelling, University of Bristol is derived from the 1989 Bangladesh Fertility Survey, [@huq:cleland:1990].
The data are a subsample of responses from 1934 women grouped in 60 districts and are available as the `contra` data set in the `MixedModels` package.

```jl
sco("contra = MixedModels.dataset(:contra)")
```

The response of interest is `use` --- whether the woman chooses to use artificial contraception.
The covariates include the district in which the woman resides, the number of live children she currently has, her age and whether she is in a rural or an urban setting.

Note that the `age` variable is centered about a particular age so some values are negative.
Regretably, the information on what the centering age was does not seem to be available.

### Plotting the binary response {#sec:plottingbinary}

Producing informative graphical displays of a binary response as it relates to covariates is somewhat more challenging that the corresponding plots for responses on a continuous scale.
If we were to plot the 1934 responses as 0/1 values versus, for example, the woman's centered age, we would end up with a rather uninformative plot because all the points would fall on one of two horizontal lines.

One approach to illustrating the structure of the data more effectively is to add *scatterplot smoother* lines
(Fig. [\[fig:Contra1\]](#fig:Contra1){reference-type="ref"
reference="fig:Contra1"})

to show the trend in the response with respect to the covariate.
Once we have the smoother lines in such a plot we can omit the data points themselves, as we did here, because they add very little information.

The first thing to notice about the plot is that the proportion of women using contraception is not linear in age, which, on reflection, makes sense.
A woman in the middle of this age range (probably corresponding to an age around 25) is more likely to use artificial contraception than is a girl in her early teens or a woman in her forties.
We also see that women in an urban setting are more likely to use contraception than those in a rural setting and that women with no live children are less likely than women who have live children.
There do not seem to be strong differences between women who have 1, 2 or 3 or more children compared to the differences between women with children and those without children.

Interestingly, the quadratic pattern with respect to age does not seem to have been noticed.
Comparisons of model fits through different software systems, as provided by the Center for Multilevel Modelling, incorporate only a linear term in age, even though the pattern is clearly nonlinear.
The lesson here is similar to what we have seen in other examples; careful plotting of the data should, whenever possible, precede attempts to fit models to the data.

### Initial GLMM fit to the contraception data {#sec:ContraGLMM}

A GLMM is fit in very much the same way as an LMM with an additional argument providing the distribution family for $\mathcal{Y}|\mathcal{B}=\mathbf{b}$.

```jl
s = """
m15 = fit(
    MixedModel,
    @formula(use ~ 1 + age + abs2(age) + livch + urban + (1|dist)),
    contra,
    Bernoulli();
    contrasts = Dict(:urban => EffectsCoding()),
    )
"""
sco(s)
```

For binary responses the distribution family to use is the [Bernoulli](https://en.wikipedia.org/wiki/Bernoulli_distribution) distribution.

Additionally a link function can be specified.
For distributions in the [exponential family](https://en.wikipedia.org/wiki/Exponential_family), the form of the distribution itself provides a *canonical link*, which is the *logit link* for the Bernoulli.

That is, the link, g, is the logit or "log-odds" function
$$
\eta=g(\mu) = \log\left(\frac{\mu}{1-\mu}\right)
$$
and the inverse link is the *logistic* function
$$
\mu=g^{-1}(\eta)=\frac{1}{1+\exp(-\eta)}
$$

The interpretation of the coefficients in this model is somewhat different from the linear mixed models coefficients that we examined previously but many of the model-building steps are similar.
A rough assessment of the utility of a particular term in the fixed-effects part of the model can be obtained from examining the estimates of the coefficients associated with it and their standard errors.
To test whether a particular term is useful we omit it from the model, refit and compare the reduced model fit to the original according to the change in deviance.

We will examine the terms in the model first and discuss the interpretation of the coefficients later in @sec:not-yet-written

Recall from Chap. that the default set of contrasts for a factor such as is offsets relative to the reference level, in this case women who do not have any live children.
Although the coefficients labeled , and are all large relative to their standard errors, they are quite close to each other.
This confirms our earlier impression that the main distinction is between women with children and those without and, for those who do have children, the number of children is not an important distinction.

After incorporating a new variable --- an indicator of whether the woman
has any children --- in the data we fit a reduced model, `m16`, with summary

```jl
s = """
contra = @transform!(DataFrame(contra), :children = :livch ≠ "0")
m16 = fit(
    MixedModel,
    @formula(use ~ 1 + age + abs2(age) + children + urban + (1|dist)),
    contra,
    Bernoulli();
    contrasts = Dict(
        :urban => EffectsCoding(),
        :children => EffectsCoding(),
    ),
)           
"""
sco(s)
```

Comparing this model to the previous model

```jl
sco("MixedModels.likelihoodratiotest(m16, m15)")
```

indicates that the reduced model is adequate.

A plot of the smoothed observed proportions versus centered age according to and
(Fig. [\[fig:Contra2\]](#fig:Contra2){reference-type="ref"
reference="fig:Contra2"})

indicates that all four groups have a quadratic trend with respect to age but the location of the peak proportion is shifted for those without children relative to those with children.
Incorporating an interaction of `age` and `children` allows for such a shift.

```jl
s = """
m17 = fit(
    MixedModel,
    @formula(use ~ 1 + age * children + abs2(age) + urban + (1|dist)),
    contra,
    Bernoulli();
    contrasts = Dict(
        :urban => EffectsCoding(),
        :children => EffectsCoding(),
    ),
)           
"""
sco(s)
```

Comparing this fitted model to the previous one

```jl
sco("MixedModels.likelihoodratiotest(m16, m17)")
```

confirms the usefulness of this term.

Continuing with the model-building we turn our attention to the random effects specification to see whether urban/rural differences vary significantly between districts and whether the distinction between childless women and women with children varies between districts.
We fit a succession of models, described in the exercises for this chapter,
before settling on model ,

```jl
s = """
m18 = fit(
    MixedModel,
    @formula(use ~ 1 + age * children + abs2(age) + urban + (1|dist&urban)),
    contra,
    Bernoulli();
    contrasts = Dict(
        :urban => EffectsCoding(),
        :children => EffectsCoding(),
    ),
) 
"""
sco(s)
```
Notice that although there are 60 distinct districts there are only 102 distinct combinations of represented in the data.
In 15 of the 60 districts there are no rural women in the sample and in 3 districts there are no urban women in the sample, as shown in


## Link functions and interpreting coefficients {#sec:GLMMlink}

To this point the only difference we have encountered between and as model-fitting functions is the need to specify the distribution family in a call to `fit`.
The formula specification is identical and the assessment of the significance of terms using likelihood ratio tests is similar.
This is intentional.
We have emphasized the use of likelihood ratio tests on terms, whether fixed-effects or random-effects terms, exactly so the approach will be general.

However, the interpretation of the coefficient estimates in the different types of models is different.
In a linear mixed model the linear predictor is the conditional mean (or "expected value") of the response given the random effects.
That is, if we assume that we know the values of the fixed-effects parameters and the random effects then the expected response for a particular combination of covariate values is the linear predictor.
Individual coefficients can be interpreted as slopes of the fitted response with respect to a numeric covariate or as shifts between levels of a categorical covariate.

To interpret the estimates of coefficients in a GLMM we must define and examine the link function that we mentioned earlier.

### The logit link function for binary responses {#sec:logitlink}

The probability model for a binary response is the Bernoulli distribution, which is about the simplest probability distribution we can concoct.
There are only two possible values: 0 and 1.
If the probability of the response 1 is $p$ then the probability of 0 must be $1-p$.
It is easy to establish that the expected value is also $p$.
For consistency across distribution families we write this expected response as $\mu$ instead of $p$.
We should, however, keep in mind that, for this distribution, $\mu$ corresponds to a probability and hence must satisfy $0\le\mu\le 1$.

In general we don't want to have restrictions on the values of the linear predictor so we equate the linear predictor to a function of $\mu$ that has an unrestricted range.
In the case of the Bernoulli distribution with the canonical link function we equate the linear predictor to the *log odds* or *logit* of the positive response.
That is 
$$\label{eq:logodds}
  \eta = \logit(\mu) = \log\left(\frac{\mu}{1-\mu}\right) .
$$

To understand why this is called the "log odds" recall that $\mu$ corresponds to a probability in $[0,1]$.
The corresponding odds ratio, $\frac{\mu}{1-\mu}$, is in $[0,\infty)$ and the logarithm of the odds ratio, $\logit(\mu)$, is in $(-\infty, \infty)$.

The inverse of the logit link function,
$$\label{eq:logitinv}
  \mu = \frac{1}{1+\exp(-\eta)}
$$
shown in

Fig. [\[fig:logitinv\]](#fig:logitinv){reference-type="ref"
reference="fig:logitinv"}, takes a value on the unrestricted range,
$(-\infty,\infty)$, and maps it to the probability range, $[0,1]$. It
happens this function is also the cumulative distribution function for
the standard logistic distribution, available in as the function . In
some presentations the relationship between the logit link and the
logistic distribution is emphasized but that often leads to questions of
why we should focus on the logistic distribution. Also, it is not clear
how this approach would generalize to other distributions such as the
Poisson or the Gamma distributions.

### Canonical link functions {#sec:canonicallink}

A way of deriving the logit link that does generalize to a class of
common distributions in what is called the *exponential family* is to
consider the logarithm of the probability function (for discrete
distributions) or the probability density function (for continuous
distributions). The probability function for the Bernoulli distribution
is $\mu$ for $y=1$ and $1-\mu$ for $y=0$. If we write this in a somewhat
peculiar way as $\mu^y+(1-\mu)^{1-y}$ for $y\in\{0,1\}$ then the
logarithm of the probability function becomes $$\label{eq:BernoulliProb}
  \log\left(\mu^y+(1-\mu)^{1-y}\right) = \log(1-\mu) +
  y\,\log\left(\frac{\mu}{1-\mu}\right) .$$ Notice that the logit link
function is the multiple of $y$ in the last term.

For members of the exponential family the logarithm of the probability
or probability density function can be expressed as a sum of up to three
terms: one that involves $y$ only, one that involves the parameters only
and the product of $y$ and a function of the parameters. This function
is the canonical link.

In the case of the Poisson distribution the probability function is
$\frac{e^{-\mu}\mu^y}{y!}$ for $y\in\{0,1,2,\dots\}$ so the log
probability function is $$\label{eq:PoissonProb}
  -\log(y!)-\mu+y\log(\mu) .$$ and the canonical link function is
$\log(\mu)$.

### Interpreting coefficient estimates {#sec:GLMMcoefficients}

Returning to the interpretation of the estimated coefficients in model
we apply exactly the same interpretation as for a linear mixed model but
taking into account that slopes or differences in levels are with
respect to the logit or log-odds function. If we wish to express results
in the probability scale then we should apply the function to whatever
combination of coefficients is of interest to us.

For example, we see from
Fig. [\[fig:Contra2\]](#fig:Contra2){reference-type="ref"
reference="fig:Contra2"} that the observed proportion of childless women
with a centered age of 0 living in a rural setting who use artificial
contraception is about 20%. The fitted value of the log-odds for a
typical district (i.e. with a random effect of zero) is corresponding to
a fitted probability of

or %.

Similarly the predicted log-odds of a childless woman with a centered
age of 0 in an urban setting of a typical district using artificial
contraception is

corresponding to a probability of

The predicted log-odds and predicted probability for a woman with
children and at the same age and location are

We should also be aware that the random effects are defined on the
linear predictor scale and not on the probability scale. A normal
probability plot of the conditional modes of the random effects for
model (Fig. [\[fig:fm13predqq\]](#fig:fm13predqq){reference-type="ref"
reference="fig:fm13predqq"})

shows that the smallest random effects are approximately -1 and the
largest are approximately 1. The numerical values and the identifier of
the combination of and for these extreme values can be obtained as

and

The exponential of the random effect is the relative odds of a woman in
a particular urban/district combination using artificial birth control
compared to her counterpart (same age, same with/without children
status, same urban/rural status) in a typical district. The odds of a
rural woman in district 1 (i.e. the value of the interaction) using
artifical contraception is

or about 40% of that of her rural counterpart in a typical district.

Notice that there is considerable variability in the lengths of the
prediction intervals in
Fig. [\[fig:fm13predqq\]](#fig:fm13predqq){reference-type="ref"
reference="fig:fm13predqq"}, unlike those from previous model fits (e.g.
Fig. [\[fig:fm01preddot\]](#fig:fm01preddot){reference-type="ref"
reference="fig:fm01preddot"}, p.  or
Fig. [\[fig:fm03ranef\]](#fig:fm03ranef){reference-type="ref"
reference="fig:fm03ranef"},
p. [\[fig:fm03ranef\]](#fig:fm03ranef){reference-type="ref"
reference="fig:fm03ranef"}). This is to be expected with data from a
highly unbalanced observational study.

Consider the cross-tabulation of counts of interviewees by district and
urban/rural status presented at the end of . The data contains responses
from 54 rural women in district 1 but only 21 rural women from district
11. Thus the bottom line in
Fig. [\[fig:fm13predqq\]](#fig:fm13predqq){reference-type="ref"
reference="fig:fm13predqq"}, from the level of the interaction, and
based on 54 responses, is shorter than the line second from the bottom,
for and based on 21 women only.
