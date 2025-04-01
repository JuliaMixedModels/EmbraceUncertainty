# Wilkinson-Rogers (1973) notation for models of (co)variance {#sec:wilkinson}

## General rules

- "Addition" (`+`) indicates additive, i.e., main effects: `a + b` indicates main effects of `a` and `b`.
- "Multiplication" (`*`) indicates crossing: main effects and interactions between two terms: `a * b` indicates main effects of `a` and `b` as well as their interaction.
- Usual algebraic rules apply (associativity and distributivity):
  - `(a + b) * c` is equivalent to `a * c + b * c`
  - `a * b * c` corresponds to main effects of `a`, `b`, and `c`, as well as all three two-way interactions and the three-way interaction.
- Categorical terms are expanded into the associated indicators/contrast variables.
- Tilde (`~`) is used to separate response from predictors.
- The intercept is indicated by `1`.
- `y ~ 1 + (a + b) * c` is read as:
  - The response variable is `y`.
  - The model contains an intercept.
  - The model contains main effects of `a`, `b`, and `c`.
  - The model contains interactions between a and c and between b and c but not a and b
- We extend this notation for mixed-effects models with the grouping notation (`|`):
  - `(1 + a | subject)` indicates "by-subject random effects for the intercept and main effect `a`".
  - This is in line with the usual statistical reading of `|` as "conditional on".

## Mixed models in Wilkinson-Rogers and mathematical notation

Models fit with MixedModels.jl are generally linear mixed-effects models with unconstrained random effects covariance matrices and homoskedastic, normally distributed residuals.
Under these assumptions, the model specification

`response ~ 1 + (age + sex) * education * n_children  + (1 | subject)`

corresponds to the statistical model

\begin{align*}
\left(Y |\mathcal{B}=b\right) &\sim N\left(X\beta + Zb, \sigma^2 I \right) \\
\mathcal{B} &\sim N\left(0, G\right)
\end{align*}

for which we wish to obtain the maximum-likelihood estimates for $G$ and thus the fixed-effects $\beta$.

- The model contains no restrictions on $G$, except that it is positive semidefinite.
- The response Y is the value of a given response.
- The fixed-effects design matrix X consists of columns for
  - the intercept, age, sex, education, and number of children (contrast coded as appropriate)
  - the interaction of all lower order terms, excluding interactions between age and sex
- The random-effects design matrix Z includes a column for
  - the intercept for each subject
