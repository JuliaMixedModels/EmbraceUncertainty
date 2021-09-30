# Overview: Conceptual Introduction to Linear Mixed-Effects Models

Linear mixed models are suitable for the analysis of multidimensional experimental and observational studies. 
In this vignette we describe the general framework and the theory behind a type of statistical model called *mixed-effects models*. 
These models are used in many different disciplines. 
Because the descriptions of the models can vary markedly between disciplines, we begin by describing what mixed-effects models are in very general terms ("The big picture of LMMs") for one type of mixed model, the *linear mixed model*

Then we illustrate for a very simple example the practice of fitting and analyzing such models using the [`MixedModels`](https://github.com/JuliaStats/MixedModels.jl) package for [`Julia`](https://julialang.org). 
This simple example allows us to illustrate the use of the `MixedModels` package for fitting such models and other functions for analyzing the fitted model. 
We also describe methods of assessing the precision of the parameter estimates and of visualizing the conditional distribution of the random effects, given the observed data.

# The big picture of linear mixed models (LMM)

## Response, covariates, and factors

LMMs, like many other types of statistical models, describe a relationship between a *response* variable and some of the *covariates* that have been measured or observed along with the response. 
The statistical model assumes that the residuals of the fitted response (i.e., not the responses) are normally -- also identically and independently -- distributed. 
This is the *first assumption* of normality in the linear mixed model. 
It is standard practice that model residuals are inspected and, if serious skew is indicated, that the response is Box-Cox transformed to fulfill this model assumption. 

In the following we distinguish between *categorical covariates* and *numerical covariates*.
Categorical covariates are  *factors*.
The important characteristic of a factor is that, for each observed value of the response, the factor takes on the value of one of a set of discrete levels.
The levels can be unordered (nominal) or ordered (ordinal). 
We use the term *covariate* when we refer to *numerical covariates*, that is to continuous measures with some distribution. In principle, statistical models are not constrained by the distribution of observations across levels of factors and covariates, but the distributions may lead to problems of model identification and they have implications for the statistical power. 

Statistical power, especially for the detection of interactions, is best for uniformly distributed factors and covariates.
In experimental designs, uniformity is achieved by balanced assignment of subjects (or other carriers of responses) to the levels of factors or combinations of factor levels.
In observational contexts, we achieve uniform distributions by stratification (e..g., on age, gender, or IQ scores).
Statistical power is also worse for skewed than normal distributions.
Therefore, although it is *not* required to meet an assumption of the statistical model, it may be useful to consider Box-Cox transformations of covariates.

## Nested and crossed random (grouping) factors

In LMMs the levels of at least one of the factors represents *units* in the data set that are assumed to be sampled, ideally randomly, from a population that is normally distributed with respect to the response. 
*This is the second assumption of normal distribution in LMMs.*  
In psychology and linguistics the observational units are often the subjects or items (e..g., texts, sentences, words, pictures) in the study.
In sociology we may look at aggregations of individuals such as we find in countries, communes, organizations, companies, schools, or families.
In agriculture the experimental units may be the plots of land or the specific plants being studied. 
We may use numbers, such as subject identifiers, to designate the particular levels that we observed but these numbers are simply labels.

Random sampling is the basis of generalization from the sample to the population.
The core statistics we will estimate in this context are variances and correlations of grand means and (quasi-)experimental effects.
These terms will be explained below. What we want to stress here is that the estimation of variances and correlations requires a larger number of units (levels) than the estimation of means.
Therefore, from a practical perspective, it is important that random factors are represented with many units.
For example, after seeing five persons on a previously unknown island will hardly allow us to get obtain reliable estimate of the variance in population height.
Large samples of subjects and words will yield reliable estimates; separate analysis for subgroups of subjects or words will substantially reduce their reliability.

When there is more than one random factor, it is useful to be clear about their relation.
The two prototypical cases are that they are *nested* or *crossed*.  In multilevel models, a special case of mixed models, the levels of the random factors are strictly nested.
For examples, at a given time, students attend classes in schools. Students, classes, and schools could be three random factors.
As soon as we look at this scenario across several school years, the nesting quickly falls apart because students may move between classes and between schools. 

In psychology and linguistics, random factors are often crossed, for example, when every subject reads every word in every sentence in a word-by-word self-paced reading experiment (or alternatively: when every word in every sentence elicits a response from every subject).
However, in an eye-movement experiment (for example), the perfect crossing on a measure like fixation duration is not attainable because of blinks or skipping of words.

In summary, the typical situation in experimental and observational studies with more than one random factor is partial crossing or partial nesting of levels of the random factors.
Linear mixed models handle these situations very well. 

## Experimental and quasi-experimental fixed factors / covariates

*Fixed experimental factor or covariate*.
In experiments the units (or levels) of the random factor(s) are assigned to manipulations implemented in their design.
The researcher controls the assignment of units of the random factor(s) (e.g., subjects, items) to experimental manipulations.
These manipulations are represented as factors with a fixed and discrete set of levels (e.g., training vs. control group) or as covariates associated with continuous numeric values (e.g., presentation times). 

*Fixed quasi-experimental factor or covariate*.
In observational studies (which can also be experiments) the units (or levels) of random factors may "bring along" characteristics that represent the levels of quasi-experimental factors or covariates beyond the control of the researcher.
Whether a a subject is female, male, or diverse or whether a word is a noun, a verb, or an adjective are examples of quasi-experimental factors of gender or word type, respectively.
Subject-related covariates are body height, body mass, and IQ scores; word-related covariates are their lengths, frequency, and cloze predictability. 

## Between-unit and within-unit factors / covariates

The distinction between between-unit and within-unit factors is always relative to a random (grouping) factor of an experimental design.
A between-unit factor / covariate is a factor for which every unit of the random factor is assigned to or characterized by only one level of the factor.
A within-unit factor is a factor for which every unit of the random factor appears at every level of the factor. 

For the typical random factor, say *Subject*, there is little ambiguity because we are used to the between-within distinction from ANOVAs, more specifically the F1-ANOVA.
In psycholinguistics, there is the tradition to test effects also for the second random factor *Item* in the F2-ANOVA.
Importantly, for a given fixed factor all four combinations are possible.
For example, *Gender* is a fixed quasi-experimental between-subject / within-item factor; word frequency is fixed quasi-experimental within-subject / between-item covariate; *Pime-target relation* is a fixed experimental  within-subject / within-item factor (assuming that targets are presented both in a primed and in an unprimed situation); and when a training manipulation is defined by the items used in the training, then in a training-control group design, the fixed factor *Group* is a fixed experimental between-subject / between-item factor.    

These distinctions are critical for setting up LMMs because variance components for (quasi-)experimental effects can only be specified for within-unit effects.
Note also that loss of data (within limits), counterbalancing or blocking of items are irrelevant for these definitions. 

## Factor-based contrasts (indicator variables) and covariate-based trends

The simplest fixed factor has two levels and the model estimates the difference between them.
When we move to factors with *k*  levels, we must decide on how we *spend* the *k-1* degrees of freedom, that is we must specify a set of contrasts.
(If we don't do it, the program chooses dummy contrasts for us.) 
The choice of contrasts determines the design or model matrix; it represents the translation of factors to contrasts (and their interactions) as indicator variables. 

The simplest specification of a covariate is to include its linear trend, that is its slope.
The slope (like a contrast) represents a difference score, that is the change in response to a one-unit change on the covariate.
For covariates we must decide on the order of the trend we want to model.

## Contrast- and trend-based fixed-effect model parameters 

Fixed factors and covariates are expected to have effects on the response. 
Fixed-effect model parameters estimate the hypothesized main and interaction effects of the study. 
The estimates of factors are based on contrasts; the estimates of covariates are based on trends.
Conceptually, they correspond to unstandardized regression coefficients in multiple regression. 

The intercept is a special regression coefficient; it estimates the value of the dependent variable when all fixed effects associated with factors and trends associated with covariates are zero.
As a rule of thumb, there is an advantage of specifying the LMM in such a way that the intercept estimates the grand mean (GM).
This happens if (a) contrasts for factors are chosen such that the intercept estimates the GM (positive: Sum, SeqDifference, or Helmert contrasts; negative: Dummy contrasts), (b) orthogonal polynomial trends are used (Helmert, anova-based), and (c) covariates are centered on their mean before inclusion in the model.
As always, there may be good theoretical reasons to depart from the default recommendation. 

The specification of contrasts / trends does not depend on the status of the fixed factor / covariate.
It does not matter whether a factor varies between or within the units of a random factor or whether it is an experimental or quasi-experimental factor.
Contrasts are *not* specified for random (grouping) factors.

## Variance components (VCs) and correlation parameters (CPs)

Variance components (VCs) and correlation parameters (CPs) are within-group model parameters; they correspond to (some of the) *within-unit* (quasi-)experimental fixed-effect model parameters.
Thus, we may be able to estimate a subject-related VC for word frequency.
If we included a linear trend for word frequency, the VC estimates the between-subject variance in these slopes.
We cannot estimate an item-related VC for the word-frequency slopes because there is only one frequency associated with words.
Analogously, we may able to estimate an item-related VC for the effect of `Gender`, but we cannot estimate a subject-related VC for this effect. 

The within-between characteristics of fixed factors and covariates relative to the random factor(s) are features of the design of the experiment or observational study.
They fundamentally constrain the specification of the LMM. 
That's why it is of upmost importance to be absolutely clear about their status.  

## Random effects (conditional means)

In this outline of the dimensions underlying the specification of an LMM, we have said nothing so far about random effects.
They figure prominently in the next sections, where we develop the mathematical foundation of LMMs.
Here we conclude with describing them as predictions for the units (levels) of the random factor(s) given the unit's data (e.g., a specific subject's data) and the model parameters. 

## Mixed-effects models

Now we turn to the mathematical foundation of *mixed-effects models* or, more simply, *mixed models*.
They are statistical models that incorporate both fixed-effects and (co-)variance parameters.
Because of the way that we will define variance parameters a model with variance parameters always includes at least one fixed-effect parameter.
Thus, any model with a variance parameter is a mixed model.

Mathematically, the random effects come into play from the outset.
We characterize the statistical model in terms of two random variables: a $q$-dimensional vector of random effects represented by the random variable $\mathcal{B}$ and an $n$-dimensional response vector represented by the random variable $\mathcal{Y}$. 
(We use upper-case “script” characters to denote random variables. 
The corresponding lower-case upright letter denotes a particular value of the random variable.) 
We observe the value, $\bf{y}$, of $\mathcal{Y}$. We do not observe the value, $\bf{b}$, of $\mathcal{B}$.

When formulating the model we describe the unconditional distribution of $\mathcal{B}$ and the conditional distribution, $(\mathcal{Y}|\mathcal{B}=\bf{b})$. 
The descriptions of the distributions involve the form of the distribution and the values of certain parameters. 
We use the observed values of the response and the covariates to estimate these parameters and to make inferences about them.

That’s the big picture.
Now let’s make this more concrete by describing a particular, versatile class of mixed models called *linear mixed models* and by studying a simple example of such a model. First we describe the data in the example.