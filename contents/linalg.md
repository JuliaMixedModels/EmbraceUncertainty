# Linear Algebra for Linear Models {#sec.LinAlg}

In this appendix we describe properties of the multivariate Gaussian (or "normal") distribution and how linear models and linear mixed models can be formulated in terms of this distribution.

We also describe some methods in numerical linear algebra that are particularly useful in working with linear models.
One of the strengths of the Julia language is the [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) package in the standard library.
The implementation of a multi-dimensional *array*, including one-dimensional *vectors* and two-dimensional *matrices*, is part of the base language.
The added value of the *LinearAlgebra* package is compact representations of special types of matrices and methods for and with matrix decompositions or factorizations.

The purpose of these descriptions is to motivate a representation of a linear mixed model that allows for fast and stable estimation of the parameters.
The estimation process requires iterative optimization of some of the parameters in the model to minimize an objective function.
Often this optimization requires hundreds or thousands of evaluations of the objective at different values of the parameters and the portion of time spent in these evaluations dominates the overall estimation time.
Thus, a fast, efficient method for evaluating the objective is crucial to making the whole process fast.

## Matrix-vector representation of linear models

A linear statistical model is often written in terms of each element of the $n$-dimensional response vector, $\mathbf{y}$, as, e.g.
$$
y_i = \beta_1 x_{i,1} + \beta_2 x_{i,2} + \dots + \beta_p x_{i,p} + \epsilon_i, \quad i=1,\dots, n
$$ {#eq:elementlinmod}
and some additional description like "where the $\epsilon_i,i=1,\dots,n$ are independently and identically distributed as $\mathcal{N}(0, \sigma^2)$".

An alternative is to write the model in terms of the $n$-dimensional *response vector*, $\mathbf{y}$, an $n\times p$ *model matrix*, $\mathbf{X}$, and a $p$-dimensional *coefficient vector*, $\boldsymbol{\beta}$, as
$$
\mathcal{Y} \sim \mathcal{N}\left(\mathbf{X}\boldsymbol{\beta}, \sigma^2\mathbf{I}\right),
$$ {#eq:mvnlinmod}
where $\mathcal{N}$ denotes the multivariate Gaussian distribution with mean $\boldsymbol{\mu}=\mathbf{X}\boldsymbol{\beta}$ and variance-covariance matrix $\boldsymbol{\Sigma}=\sigma^2\mathbf{I}$.
(In what follows we will refer to the *variance-covariance matrix* as simply the *covariance matrix*.)

Before considering properties of and computational methods for the model @eq:mvnlinmod we will describe some of the properties of the [multivariate Gaussian distribution](https://en.wikipedia.org/wiki/Multivariate_normal_distribution).


## The multivariate Gaussian distribution

Just as a univariate Gaussian distribution can be written by specifying the (scalar) mean, $\mu$, and the variance, $\sigma^2$, as $\mathcal{N}(\mu, \sigma^2)$, a multivariate Gaussian distribution is characterized by its $n$-dimensional mean vector, $\boldsymbol{\mu}$, and its $n\times n$ variance-covariance matrix, $\boldsymbol{\Sigma}$, as $\mathcal{N}(\boldsymbol{\mu}, \boldsymbol{\Sigma})$.

The density function for a univariate Gaussian distribution is the familiar "bell curve"
$$
f(x; \mu, \sigma^2)=\frac{1}{\sqrt{2\pi\sigma^2}}\exp\left(\frac{-\left(x-\mu\right)^2}{2\sigma^2}\right)
$$ {#eq:univariateGaussianpdf}
and probabilities defined by this density are most easily evaluated by standardizing the deviation, $x-\mu$, as $z=\frac{x-\mu}{\sigma}$.
(This is why $\sigma$ is called the *standard deviation*.)

To be able to evaluate $\sigma$, the variance, $\sigma^2$, must be positive, or at least non-negative.
If $\sigma^2=0$ then all the probability is concentrated at a single point, $x=\mu$, and we no longer have a probability density, in the usual way of thinking of one.
The density shrinks to a point mass and the distribution is said to be [degenerate](https://en.wikipedia.org/wiki/Degenerate_distribution).

Similar constraints apply to the covariance matrix, $\boldsymbol{\Sigma}$.
Because the covariance of the i'th and j'th elements does not depend upon the order in which we write them, $\boldsymbol\Sigma$ must be symmetric.
That is,
$$
\boldsymbol{\Sigma}' = \boldsymbol{\Sigma}
$$ {#eq:Sigmasym}
Furthermore, to define a proper multivariate density, $\boldsymbol\Sigma$ must be [positive definite](https://en.wikipedia.org/wiki/Definite_matrix), which means that for any non-zero vector, $\mathbf{x}$, the *quadratic form* defined by $\boldsymbol\Sigma$ must be positive.
That is
$$
  \mathbf{x}'\boldsymbol{\Sigma}\mathbf{x}>0,\quad\forall\,\mathbf{x}\ne\mathbf{0} .
$$ {#eq:positiveDef}
(the symbol $\forall$ means "for all").
Positive definiteness implies that the *precision matrix*, $\boldsymbol\Sigma^{-1}$, exists and is also positive definite.
It also implies that there are "matrix square roots" of $\boldsymbol\Sigma$ in the sense that there are matrices $\mathbf{A}$ such that $\mathbf{A}'\mathbf{A}=\boldsymbol{\Sigma}$.
(The reason for writing $\mathbf{A}'\mathbf{A}$ and not simply the square of $\mathbf{A}$ is that $\mathbf{A}$ is not required to be symmetric but $\mathbf{A}'\mathbf{A}$ will be symmetric, even in $\mathbf{A}$ is not.)

One such "square root" of a positive definite $\boldsymbol\Sigma$ is the [Cholesky factor](https://en.wikipedia.org/wiki/Cholesky_decomposition), which corresponds to $n\times n$ upper-triangular matrix, $\mathbf{R}$, such that
$$
\boldsymbol{\Sigma}=\mathbf{R}'\mathbf{R} .
$$ {#eq:upperCholesky}
This factor is usually called $\mathbf{R}$ because it appears without the transpose as the right-hand multiplicant in @eq:upperCholesky.
An alternative expression is written with the lower-triangular $\mathbf{L}$ on the left as
$$
\boldsymbol{\Sigma}=\mathbf{L}\mathbf{L}',
$$ {#eq:lowerCholesky}
with the obvious relationship that $\mathbf{L}=\mathbf{R}'$.
To add to the confusion, the `cholesky` function in the *LinearAlgebra* package produces a factorization where the lower-triangular factor on the left is called `L` and the upper-triangular factor on the right is called `U`.

The factor $\mathbf{R}$ or $\mathbf{L}$ can be evaluated directly from the elements of $\boldsymbol\Sigma$.
For example, the non-zeros in the first two rows of $\mathbf{L}$ are evaluated as
$$
\begin{align}
  \mathbf{L}_{1,1}&=\sqrt{\boldsymbol{\Sigma}_{1,1}}\\
  \mathbf{L}_{2,1}&=\boldsymbol{\Sigma}_{2,1}/\mathbf{L}_{1,1}\\
  \mathbf{L}_{2,2}&=\sqrt{\boldsymbol{\Sigma}_{2,2}-\mathbf{L}_{2,1}^2}
\end{align}
$$ {#eq:lowerCholtworows}
Evaluating the diagonal elements involves taking a square root.
By convention we choose the positive square root for the Cholesky factor with the result that the diagonal elements of $\mathbf{L}$ are all positive.

### Some properties of triangular matrices

A triangular matrix with non-zero diagonal elements is non-singular.
One way to show this is because its [determinant](https://en.wikipedia.org/wiki/Determinant), written $\left|\mathbf{L}\right|$, which is the product of its diagonal elements, is non-zero.
In the case of a Cholesky factor the determinant will be positive because all the diagonal elements are positive.

A more straightforward way of showing that such a matrix is non-singular is to show how a triangular system of equations, like
$$
\mathbf{Lx}=\mathbf{b}
$$ {#eq:triangularsystem}
can be solved.
In the case of a lower-triangular system the method is called *forward solution*, with the sequence of scalar equations
$$
\begin{align}
x_1&=b_1/\mathbf{L}_{1,1}\\
x_2&=\left(b_2-x_1\mathbf{L}_{2,1}\right)/\mathbf{L}_{2,2}\\
x_3&=\left(b_3-x_1\mathbf{L}_{3,1}-x_2\mathbf{L}_{3,2}\right)/\mathbf{L}_{3,3}
\end{align}
$$ {#eq:forwardsolve}
and so on.

One point to note here is that $b_1$ is not needed after $x_1$ is evaluated, $b_2$ is not needed after $x_2$ is evaluated, and so on.
That is, the forward solution can be carried out *in place* with each element of $\mathbf{b}$ overwriting the corresponding element of $\mathbf{x}$.
This property is useful for avoiding allocation of storage in each evaluation of the objective function.

The corresponding method of solving an upper-triangular system of equations is called *backward solution*, where $b_n$ is evaluated first, then $b_{n-1}$, and so on.

Repeated forward solution (or backward solution for upper triangular) can be used to evaluate the inverse, $\mathbf{L}^{-1}$, of a lower triangular matrix, $\mathbf{L}$.
However, a general rule in numerical linear algebra is that you rarely need to evaluate the full inverse of a matrix.
Solving a triangular system like @eq:triangularsystem by evaluating $\mathbf{L}^{-1}$ and forming the product 
$$
\mathbf{x} = \mathbf{L}^{-1}\mathbf{b}
$$ {#eq:inversesolve}
involves doing roughly $n$ times as much work as solving the system directly, as in @eq:forwardsolve.
Requiring that the inverse of a matrix must be evaluated to solve a linear system is like saying that a quotient, $a/b$, must be evalated by calculating $b^{-1}$, the reciprocal of $b$, then evaluating the product $b^{-1}a$, instead of evaluating the quotient directly.

In a derivation we may write an expression like $\mathbf{L}^{-1}\mathbf{b}$ but the evaluation is performed by solving a system like @eq:forwardsolve.

### Positive definiteness and the Cholesky factor

It turns out that the ability to form the Cholesky factor, which means that all the quantities like $\boldsymbol{\Sigma}_{2,2}-\mathbf{L}_{2,1}^2$, whose square roots form the diagonal of $\mathbf{L}$, evaluate to positive numbers, is equivalent to $\boldsymbol\Sigma$ being positive definite.
It is straightforward to show that having a Cholesky factor implies that $\boldsymbol\Sigma$ is positive definite, because
$$
\mathbf{x}'\boldsymbol{\Sigma}\mathbf{x} = \mathbf{x}'\mathbf{R}'\mathbf{R}\mathbf{x}=\left(\mathbf{Rx}\right)'\mathbf{Rx}=\left\|\mathbf{Rx}\right\|^2
$$ {#eq:cholimpliesposdef}
where $\left\|\mathbf{v}\right\|^2$ is the squared length of the vector $\mathbf{v}$.
Because $\mathbf{R}$ is non-singular, $\mathbf{x}\ne\mathbf{0}\implies\mathbf{Rx}\ne\mathbf{0}$ and the squared length in @eq:cholimpliesposdef is greater than zero.

The other direction is a bit more complicated to prove but essentially it amounts to showing that if the process of the generating the Cholesky factor requires the square root of a non-positive number to obtain a diagonal element then there is a direction in which the quadratic form gives a non-positive result.

In practice, the easiest way to check a symmetric matrix to see if it is positive definite is to attempt to evaluate the Cholesky factor and check whether that succeeds.
This is exactly what the `isposdef` methods in the *LinearAlgebra* package do.

### Density of the multivariate Gaussian

For the general multivariate normal distribution, $\mathcal{N}(\boldsymbol{\mu},\boldsymbol{\Sigma})$, where $\boldsymbol{\Sigma}$ is positive definite with lower Cholesky factor $\mathbf{L}$, the probability density function is
$$
\begin{align}
f(\mathbf{x};\boldsymbol{\mu},\boldsymbol{\Sigma})&=
\frac{1}{\sqrt{(2\pi)^n\left|\boldsymbol{\Sigma}\right|}}
\exp\left(\frac{-[\mathbf{x}-\boldsymbol{\mu}]'\boldsymbol{\Sigma}^{-1}[\mathbf{x}-\boldsymbol{\mu}]}{2}\right)\\
&=\frac{1}{\sqrt{(2\pi)^n}\left|\mathbf{L}\right|}
\exp\left(\frac{-[\mathbf{x}-\boldsymbol{\mu}]'{\mathbf{L}'}^{-1}\mathbf{L}^{-1}[\mathbf{x}-\boldsymbol{\mu}]}{2}\right)
\end{align}
$$ {#eq:mvndensity}
and the standardizing transformation becomes
$$ 
\mathbf{z}=\mathbf{L}^{-1}[\mathbf{x}-\boldsymbol{\mu}] .
$$ {#eq:mvnstandardizing}
which, in practice, means using forward solution on the lower-triangular system of equations
$$
\mathbf{Lz}=\mathbf{x}-\boldsymbol{\mu}
$$ {#eq:mvnstandardizingsol}

Note that the standardizing transformation gives us a way to simulate values from a general $n$-dimensional multivariate Gaussian, $\mathcal{X}\sim\mathcal{N}(\boldsymbol{\mu},\boldsymbol{\Sigma})$ as
$$
\mathbf{x}=\boldsymbol{\mu}+\mathbf{L}\mathbf{z}
$$ {#eq:simulatemvn}

where $\mathbf{z}$ is simulated from the $n$-dimensional *standard multivariate Gaussian*, $\mathcal{Z}\sim\mathcal{N}(\mathbf{0},\mathbf{I})$, which is $n$ independent univariate standard normal distributions.

### Linear functions of a multivariate Gaussian

In general, if $\mathcal{X}$ is an $n$-dimensional random variable with mean $\boldsymbol{\mu}$ and covariance matrix $\boldsymbol{\Sigma}$ and $\mathbf{A}$ is a matrix with $n$ columns then the mean and variance of $\mathcal{U}=\mathbf{A}\mathcal{X}$ are given by
$$
\require{unicode}
𝔼\left[\mathcal{U}\right] =
𝔼\left[\mathbf{A}\mathcal{X}\right] =
\mathbf{A}𝔼\left[\mathcal{X}\right] =
\mathbf{A}\boldsymbol{\mu}
$$ {#eq:mvexpectedlin}
and
$$
\begin{align}
\text{Var}\left(\mathcal{U}\right)
&=𝔼\left[\left(\mathcal{U}-𝔼\left[\mathcal{U}\right]\right)\left(\mathcal{U}-𝔼\left[\mathcal{U}\right]\right)'\right]\\
&=𝔼\left[\left(\mathbf{A}\mathcal{X}-\mathbf{A}\boldsymbol{\mu}\right)\left(\mathbf{A}\mathcal{X}-\mathbf{A}\boldsymbol{\mu}\right)'\right]\\
&=𝔼\left[\mathbf{A}\left(\mathcal{X}-\boldsymbol{\mu}\right)\left(\mathcal{X}-\boldsymbol{\mu}\right)'\mathbf{A}'\right]\\
&=\mathbf{A}\,𝔼\left[\left(\mathcal{X}-\boldsymbol{\mu}\right)\left(\mathcal{X}-\boldsymbol{\mu}\right)'\right]\mathbf{A}'\\
&=\mathbf{A}\text{Var}(\mathcal{X})\mathbf{A}'\\
&=\mathbf{A}\boldsymbol{\Sigma}\mathbf{A}'
\end{align}
$$ {#eq:mvvarlin}

A linear function, $\mathcal{U}=\mathbf{A}\mathcal{X}$, of a multivariate Gaussian distribution, $\mathcal{X}\sim\mathcal{N}(\boldsymbol{\mu},\boldsymbol{\Sigma})$, is also Gaussian and these relationships imply that
$$
\mathcal{U}\sim\mathcal{N}(\mathbf{A}\boldsymbol{\mu}, \mathbf{A}\boldsymbol{\Sigma}\mathbf{A}')
$$ {#eq:mvnlinfunc}

For the special case of $\mathbf{A}$ being of dimension $1\times n$ (i.e. a *row vector*), the expression for the $1\times 1$ covariance matrix is the quadratic form defined by $\boldsymbol{\Sigma}$, which is why $\boldsymbol{\Sigma}$ must be positive definite for the conditional distributions to be non-degenerate.

## Back at the linear model

The probability density function for the linear model, @eq:mvnlinmod, is
$$
\begin{align}
f(\mathbf{y}; \boldsymbol{\beta}, \sigma^2)&=
\frac{1}{\sqrt{2\pi\left|\sigma^2\mathbf{I}\right|}}
\exp\left(\frac{-[\mathbf{y}-\mathbf{X}\boldsymbol{\beta}]'
\left(\sigma^2\mathbf{I}\right)^{-1}[\mathbf{y}-\mathbf{X}\boldsymbol{\beta}]}{2}\right)\\
&=\left(2\pi\sigma^2\right)^{-n/2}\exp\left(-\left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2/\left(2\sigma^2\right)\right)
\end{align}
$$ {#eq:linmoddensity}

@eq:linmoddensity describes the density of the random variable, $\mathcal{Y}$, representing the observations, given the values of the parameters, $\boldsymbol{\beta}$ and $\sigma^2$.
For parameter estimation we use the *likelihood function*, which is the same expression as @eq:linmoddensity but regarded as function of the parameters, $\boldsymbol{\beta}$ and $\sigma^2$, with the observed response, $\mathbf{y}$, fixed.
$$
L(\boldsymbol{\beta},\sigma^2;\mathbf{y})=
\left(2\pi\sigma^2\right)^{-n/2}\exp\left(-\left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2/\left(2\sigma^2\right)\right)
$$ {#eq:linmodlikelihood}
The [maximum likelihood](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation) estimates of the parameters are, as the name implies, the values of $\boldsymbol{\beta}$ and $\sigma^2$ that maximize the expression on the right of @eq:linmodlikelihood .

Because the logarithm is a [monotone increasing](https://en.wikipedia.org/wiki/Monotonic_function) function, the maximum likelihood estimates will also maximize the *log-likelihood*
$$
\begin{align}
\ell(\boldsymbol{\beta},\sigma^2;\mathbf{y})
&=\log L(\boldsymbol{\beta},\sigma^2;\mathbf{y})\\
&=-\frac{n}{2}\log(2\pi\sigma^2)-\frac{\left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2}{2\sigma^2}
\end{align}
$$ {#eq:linmodloglike}
Usually the log-likelihood is easier to optimize, either algebraically or numerically, than the likelihood itself.

To avoid the negative signs and the factors of 2 in the denominator, we often convert the log-likelihood to the *deviance scale*, which is negative twice the log-likelihood,
$$
\begin{align}
d(\boldsymbol{\beta},\sigma^2;\mathbf{y})
&=-2\ell(\boldsymbol{\beta},\sigma^2; \mathbf{y})\\
&=n\log(2\pi\sigma^2)+\frac{\left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2}{\sigma^2} .
\end{align}
$$ {#eq:linmoddevscale}
Because of the negative sign, the maximum likelihood estimates are those that *minimize* $d(\boldsymbol{\beta},\sigma^2;\mathbf{y})$.

(The term *deviance scale* is used for $d(\boldsymbol{\beta},\sigma^2;\mathbf{y})$ rather than [deviance](https://en.wikipedia.org/wiki/Deviance_(statistics)) because the deviance involves an additive shift, which is a correction for the saturated model - see the link.
It is obvious what the saturated model should be for the linear model but not for the linear mixed model so, to avoid confusion, we refer to the log-likelihood on the deviance scale as the *objective*.)

The form of @eq:linmoddevscale makes it easy to determine the maximum likelihood estimates.
Because $\boldsymbol{\beta}$ appears only in the sum of squared residuals expression, $\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\|^2$, we minimize that with respect to $\boldsymbol{\beta}$
$$
\widehat{\boldsymbol{\beta}}=
\arg\min_{\boldsymbol{\beta}}\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\|^2 ,
$$ {#eq:leastsquaresest}
where $\arg\min_{\boldsymbol{\beta}}$ means the value of $\boldsymbol{\beta}$ that minimizes the expression that follows.

Let $r^2(\widehat{\boldsymbol{\beta}}) = \left\|\mathbf{y}-\mathbf{X}\widehat{\boldsymbol{\beta}}\right\|^2$ be the minimum sum of squared residuals.
Substituting this value into @eq:linmoddevscale, differentiating with respect to $\sigma^2$, and setting this derivative to zero gives
$$
\widehat{\sigma^2}=\frac{r^2(\widehat{\boldsymbol{\beta}})}{n}
$$

## Minimizing the sum of squared residuals

A condition for $\widehat{\boldsymbol{\beta}}$ to minimize the sum of squared residuals is that the *gradient*
$$
\nabla r^2(\boldsymbol{\beta})=-2\mathbf{X}'(\mathbf{y}-\mathbf{X}\boldsymbol{\beta})
$$ {#eq:sumsqgrad}
be zero at $\widehat{\boldsymbol{\beta}}$.
This condition can be rewritten as
$$
\mathbf{X}'\mathbf{X}\widehat{\boldsymbol{\beta}}=\mathbf{X}'\mathbf{y} ,
$$ {#eq:normaleq}
which are called the *normal equations*.

The term *normal* in this expression comes from the fact that requiring the gradient, @eq:sumsqgrad, to be zero is equivalent to requiring that the *residual vector*, $\mathbf{y}-\mathbf{X}\widehat{\boldsymbol{\beta}}$, be perpendicular, or *normal*, to the columns of $\mathbf{X}$.

When the model matrix, $\mathbf{X}$, is of *full column rank*, which means
$$
\mathbf{X}\boldsymbol{\beta}\ne\mathbf{0}\quad\forall\boldsymbol{\beta}\ne\mathbf{0} ,
$$ {#eq:fullcolumnrank}
then the quadratic form defined by $\mathbf{X}'\mathbf{X}$ is positive definite and has a Cholesky factor, say $\mathbf{R}_{XX}$, and the normal equations can be solved in two stages.
First, solve
$$
\mathbf{R}_{XX}'\mathbf{r}_{Xy}=\mathbf{X}'\mathbf{y}
$$ {#eq:rXydef}
for $\mathbf{r}_{Xy}$ using forward solution, then solve
$$
\mathbf{R}_{XX}\widehat{\boldsymbol{\beta}}=\mathbf{r}_{Xy}
$$ {#eq:betahatchol}
for $\widehat{\boldsymbol{\beta}}$ using backward solution.

An alternative approach is to write the residual sum of squares as a quadratic form
$$
\begin{align}
r^2(\boldsymbol{\beta})&=\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\|^2\\
&=(\mathbf{y}-\mathbf{X}\boldsymbol{\beta})'(\mathbf{y}-\mathbf{X}\boldsymbol{\beta})\\
&=(\mathbf{X}\boldsymbol{\beta}-\mathbf{y})'(\mathbf{X}\boldsymbol{\beta}-\mathbf{y})\\
&=\begin{bmatrix}\boldsymbol{\beta}&-1\end{bmatrix}
\begin{bmatrix}
  \mathbf{X}'\mathbf{X} & \mathbf{X}'\mathbf{y}\\
  \mathbf{y}'\mathbf{X} & \mathbf{y}'\mathbf{y}
\end{bmatrix}
\begin{bmatrix}
  \boldsymbol{\beta}\\
  -1
\end{bmatrix}\\
&=\begin{bmatrix}\boldsymbol{\beta}&-1\end{bmatrix}
\begin{bmatrix}
  \mathbf{R}_{XX}' & \mathbf{0}\\
  \mathbf{r}_{Xy}' & r_{yy}
\end{bmatrix}
\begin{bmatrix}
  \mathbf{R}_{XX} & \mathbf{r}_{Xy}\\
  \mathbf{0} & r_{yy}
\end{bmatrix}
\begin{bmatrix}
  \boldsymbol{\beta}\\
  -1
\end{bmatrix}\\
&=\left\|
\begin{bmatrix}
  \mathbf{R}_{XX} & \mathbf{r}_{Xy}\\
  \mathbf{0} & r_{yy}
\end{bmatrix}
\begin{bmatrix}
  \boldsymbol{\beta}\\
  -1
\end{bmatrix}\right\|^2\\
&=\left\|\mathbf{R}_{XX}\boldsymbol{\beta}-\mathbf{r}_{Xy}\right\|^2
+ r_{yy}^2
\end{align}
$$ {#eq:extendedqf}

The first term, $\left\|\mathbf{R}_{XX}\boldsymbol{\beta}-\mathbf{r}_{Xy}\right\|^2$, is non-negative and can be made zero by solving @eq:betahatchol for $\widehat{\boldsymbol{\beta}}$.
Thus, the minimum sum of squared residuals is $r_{yy}^2$.

One consequence of this derivation is that the minimum sum of squared residuals can be evaluated directly from the extended Cholesky factor
$$
\begin{bmatrix}
  \mathbf{R}_{XX} & \mathbf{r}_{Xy}\\
  \mathbf{0} & r_{yy}
\end{bmatrix}
$$ {#eq:extendedcholfac}
without needing to solve for $\widehat{\boldsymbol{\beta}}$ first.
This is not terribly important for a linear model where the evaluation of $\widehat{\boldsymbol{\beta}}$ and the residual is typically done only once.
However, for the linear mixed model, a similar calculation must be done for every evaluation of the objective in the iterative optimization, and being able to evaluate the minimum penalized sum of squared residuals without solving for parameter values and without needing to evaluate the residual saves a non-negligible amount of time and effort.

## Numerical example

Suppose we wish to fit a simple linear regression model to the reaction time as a function of days of sleep deprivation to the data from subject `S372` in the `sleepstudy` dataset.
```jl
s = """
S372 = last(groupby(DataFrame(MixedModels.dataset(:sleepstudy)), :subj))
"""
sco(s; process=without_caption_label)
```
The model matrix and the response vector can be constructed as
```jl
sco("X = hcat(ones(nrow(S372)), S372.days)")
```
and
```jl
sco("y = S372.reaction")
```
from which we obtain the Cholesky factor
```jl
sco("chfac = cholesky!(X'X)")
```
(Recall that the upper triangular Cholesky factor is the `U` property of the `Cholesky` type.)

The `\` operator with a `Cholesky` factor on the left performs both the forward and backward solutions to obtain the least squares estimates
```jl
sco("β̂ = chfac\\(X'y)")
```

Alternatively, we could carry out the two solution of triangular systems explicitly by first solving for $\mathbf{r}_{Xy}$
```jl
sco("rXy = ldiv!(chfac.L, X'y)")
```
then solving in-place to obtain $\widehat{\boldsymbol{\beta}}$
```jl
sco("ldiv!(chfac.U, rXy)")
```

The residual vector, $\mathbf{y}-\mathbf{X}\widehat{\boldsymbol{\beta}}$, is
```jl
sco("r = y - X * β̂")
```
with geometric length or "norm",
```jl
sco("norm(r)")
```

For the extended Cholesky factor, create the extended matrix of sums of squares and cross products
```jl
s = """
crprod = let x = S372.days
    Symmetric(
        [length(x) sum(x)       sum(y)
             0.    sum(abs2, x) dot(x, y)
             0.      0.         sum(abs2, y)
        ],
        :U
    )
end
"""
sco(s)
```
The call to `Symmetric` with the second argument the symbol `:U` indicates that the matrix should be treated as symmetric but only the upper triangle is given.

The Cholesky factor of the `crprod` reproduces $\mathbf{R}_{XX}$, $\mathbf{r}_{Xy}$, and the norm of the residual, $r_{yy}$.
```jl
sco("extchfac = cholesky(crprod)")
```
and information from which the parameter estimates can be evaluated.
```jl
s = """
β̂ ≈ ldiv!(
  UpperTriangular(view(extchfac.U, 1:2, 1:2)),
  copy(view(extchfac.U, 1:2, 3)),
)
"""
sco(s)
```
The operator `≈` is a check of approximate equality of floating point numbers or arrays.
Exact equality of floating point results from "equivalent" calculations cannot be relied upon.

```jl
sco("norm(r) ≈ extchfac.U[3,3]")
```

## Alternative decompositions of X

There are two other decompositions of the model matrix $\mathbf{X}$ or the augmented model matrix $[\mathbf{X,y}]$ that can be used to evaluate the least squares estimates; the [QR decomposition](https://en.wikipedia.org/wiki/QR_decomposition) and the [singular value decomposition (SVD)](https://en.wikipedia.org/wiki/Singular_value_decomposition).

The QR decomposition expresses $\mathbf{X}$ as the product of an *orthogonal* matrix, $\mathbf{Q}$, and an upper triangular matrix $\mathbf{R}$.
The upper triangular $\mathbf{R}$ is related to the upper triangular Cholesky factor in that the numerical values are the same but the signs can be different.
In particular, the usual way of creating $\mathbf{Q}$ and $\mathbf{R}$ using [Householder transformations](https://en.wikipedia.org/wiki/Householder_transformation) typically results in the first row of $\mathbf{R}$ from the `qr` function being the negative of the first row of the upper Cholesky factor.

```jl
s = """
qrfac = qr(X);
qrfac.R
"""
sco(s)
```

Just as the Cholesky factor can be used on the left of the `\` operator, so can the `qr` factor but with `y` on the right.

```jl
sco(raw"b3 = qrfac\y")
``` 

The matrix $\mathbf{R}$ is returned as a square matrix with the same number of columns as $\mathbf{X}$.
That is, if $\mathbf{X}$ is of size $n\times p$ where $n>p$, as in the example, then $\mathbf{R}$ is $p\times p$, as shown above.

The matrix $\mathbf{Q}$ is usually considered to be an $n\times n$ orthogonal matrix, which means that its transpose is its inverse
$$
\mathbf{Q'Q}=\mathbf{QQ'}=\mathbf{I}
$$ {#eq:orthogonalQ}
To form the product $\mathbf{QR}$ the matrix $\mathbf{R}$ is treated as if it were $n\times p$ with zeros below the main diagonal.

The $n\times n$ matrix $\mathbf{Q}$ can be very large if $n$, the number of observations, is large but it does not need to be explicitly evaluated.
In practice $\mathbf{Q}$ is a "virtual" matrix represented as a product of Householder reflections that only require storage of the size of $\mathbf{X}$.
The effect of multiplying a vector or matrix by $\mathbf{Q}$ or by $\mathbf{Q}'$ is achieved by applying the Householder reflections in a particular order.

```jl
sco("rXy2 = qrfac.Q'y")
```

```jl
sco("b4 = ldiv!(UpperTriangular(qrfac.R), rXy2[1:2])")
```

Forming the QR decomposition is a direct, non-iterative, calculation, like forming the Cholesky factor.
Forming the SVD, by contrast, is usually an iterative calculation.
(It should be noted that modern methods for evaluating the SVD are very fast for an iterative calculation.)
The SVD consists of two orthogonal matrices, the $n\times n$ $\mathbf{U}$ and the $p\times p$ $\mathbf{V}$ and an $n\times p$ matrix $\mathbf{S}$ that is zero off the main diagonal, where
$$
\mathbf{X}=\mathbf{USV'} .
$$

Unlike the $\mathbf{Q}$ in the QR decomposition, the orthogonal matrices $\mathbf{U}$ and $\mathbf{V}$ are explicitly evaluated.
Because of this, the default for the `svd` function is to produce a compact form where $\mathbf{U}$ is $n\times p$ and only the diagonal of $\mathbf{S}$ is returned.
```jl
sco("Xsvd = svd(X)")
```

If all the singular values are non-zero, as is the case here, the least squares solution $\widehat{\boldsymbol{\beta}}$ can be obtained as
$$
\mathbf{V}\mathbf{S}^{-1}\mathbf{U}'\mathbf{y}
$$ {#eq:pseudoinv}
for the diagonal $\mathbf{S}$.

```jl
sco("b5 = Xsvd.V * (Xsvd.U'y ./ Xsvd.S)")
```

In the extensions to linear mixed-effects models we will emphasize the Cholesky factorization over the QR decomposition or the SVD.

## Linear mixed-effects models

As described in @bates.maechler.etal:2015 , a linear mixed-effects model is based on two vector-valued random variables: the $q$-dimensional vector of random effects, $\mathcal{B}$, and the $n$-dimensional response vector, $\mathcal{Y}$. @eq:LMMdist defines the unconditional distribution of $\mathcal{B}$ and the conditional distribution of $\mathcal{Y}$, given $\mathcal{B}=\mathbf{b}$, as multivariate Gaussian distributions of the form
$$
\begin{aligned}
  (\mathcal{Y}|\mathcal{B}=\mathbf{b})&\sim\mathcal{N}(\mathbf{X}\boldsymbol{\beta}+\mathbf{Z}\mathbf{b},\sigma^2\mathbf{I})\\
  \mathcal{B}&\sim\mathcal{N}(\mathbf{0},\boldsymbol{\Sigma}_\theta) .
\end{aligned}
$$

The $q\times q$, symmetric, variance-covariance matrix, $\mathrm{Var}(\mathcal{B})=\boldsymbol{\Sigma}_\theta$, depends on the *variance-component parameter vector*, $\boldsymbol{\theta}$, through a lower triangular *relative covariance factor*, $\Lambda_\theta$ as
$$
\boldsymbol{\Sigma}_\theta=\sigma^2\boldsymbol{\Lambda}_\theta\boldsymbol{\Lambda}_\theta' .
$$
(Recall that the lower Cholesky factor is generally written $\mathbf{L}$.
In this case the lower Cholesky factor contains parameters and is named with the corresponding Greek letter, $\Lambda$.)

Many computational formulas for linear mixed models are written in terms of the *precision matrix*, $\boldsymbol{\Sigma}_\theta^{-1}$.
Such formulas will become unstable as $\boldsymbol{\Sigma}_\theta$ approaches singularity.
And it can do so.
It is a fact that singular (i.e. non-invertible) $\boldsymbol{\Sigma}_\theta$ can and do occur in practice, as we have seen in some of the examples in earlier chapters.
Moreover, during the course of the numerical optimization by which the parameter estimates are determined, it is frequently the case that the deviance or the REML criterion will need to be evaluated at values of $\boldsymbol{\theta}$ that produce a singular $\boldsymbol{\Sigma}_\theta$.
Because of this we will take care to use computational methods that can be applied even when $\boldsymbol{\Sigma}_\theta$ is singular and are stable as $\boldsymbol{\Sigma}_\theta$ approaches singularity.

As defined in @eq:relcovfac, a relative covariance factor, $\Lambda_\theta$, is any matrix that satisfies
$$
\boldsymbol{\Sigma}_\theta=\sigma^2\Lambda_\theta\Lambda_\theta' .
$$
According to this definition, $\boldsymbol{\Sigma}$ depends on both $\sigma$ and $\theta$, and we should write it as $\boldsymbol{\Sigma}_{\sigma,\theta}$.
However, we will blur that distinction and continue to write $\text{Var}(\mathcal{B})=\boldsymbol{\Sigma}_\theta$.
Another technicality is that the *common scale parameter*, $\sigma$, could, in theory, be zero.
We will show that in practice the only way for its estimate, $\widehat{\sigma}$, to be zero is for the fitted values from the fixed-effects only, $\mathbf{X}\widehat{\boldsymbol{\beta}}$, to be exactly equal to the observed data.
This occurs only with data that have been (incorrectly) simulated without error.
In practice we can safely assume that $\sigma>0$.
However, $\Lambda_\theta$, like $\boldsymbol{\Sigma}_\theta$, can be singular.

The computational methods in the *MixedModels* package are based on $\Lambda_\theta$ and do not require evaluation of $\boldsymbol{\Sigma}_\theta$.
In fact, $\boldsymbol{\Sigma}_\theta$ is explicitly evaluated only at the converged parameter estimates.

The spherical random effects, $\mathcal{U}\sim\mathcal{N}(\mathbf{0},\sigma^2\mathbf{I}_q)$, determine $\mathcal{B}$ as
$$
  \mathcal{B}=\Lambda_\theta\mathcal{U} .
$$ {#eq:sphericalRE}
Although it may seem more intuitive to write $\mathcal{U}$ as a linear transformation of $\mathcal{B}$, we cannot do that when $\Lambda_\theta$ is singular, which is why @eq:sphericalRE is in the form shown.

We can easily verify that @eq:sphericalRE provides the desired distribution for $\mathcal{B}$.
As a linear transformation of a multivariate Gaussian random variable, $\mathcal{B}$ will also be multivariate Gaussian with mean
$$
𝔼\left[\mathcal{B}\right]=
𝔼\left[\boldsymbol{\Lambda}_\theta\mathcal{U}\right]=
\boldsymbol{\Lambda}_\theta\,𝔼\left[\mathcal{U}\right]=
\boldsymbol{\Lambda}_\theta\mathbf{0}=\mathbf{0}
$$
and covariance matrix
$$
\text{Var}(\mathcal{B})=
\boldsymbol{\Lambda}_\theta\text{Var}(\mathcal{U})\boldsymbol{\Lambda}\theta'=
\sigma^2\boldsymbol{\Lambda}_\theta\boldsymbol{\Lambda}_\theta'=\boldsymbol{\Sigma}_\theta
$$

Just as we concentrate on how $\boldsymbol{\theta}$ determines $\Lambda_\theta$, not $\boldsymbol{\Sigma}_\theta$, we will concentrate on properties of $\mathcal{U}$ rather than $\mathcal{B}$.
In particular, we now define the model according to the distributions
$$
  \begin{aligned}
  (\mathcal{Y}|\mathcal{U}=\mathbf{u})&\sim\mathcal{N}(\mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\beta,\sigma^2\mathbf{I}_n)\\
  \mathcal{U}&\sim\mathcal{N}(\mathbf{0},\sigma^2\mathbf{I}_q) .
  \end{aligned}
$$ {#eq:condYgivenU}

The joint density for $\mathcal{Y}$ and $\mathcal{U}$ is the product of densities of the two distributions shown in @eq:condYgivenU.
That is
$$
f_{\mathcal{Y},\mathcal{U}}(\mathbf{y},\mathbf{u})=
\frac{1}{\left(2\pi\sigma^2\right)^{-(n+q)/2}}\exp
\left(\frac{\left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}
-\mathbf{Z}\boldsymbol{\Lambda}_\theta\mathbf{u}\right\|^2+
\left\|\mathbf{u}\right\|^2}{-2\sigma^2}\right) .
$$ {#eq:YUjointdensity}

To evaluate the likelihood for the parameters, $\boldsymbol{\theta}$, $\boldsymbol{\beta}$, and $\sigma^2$, given the observed response, $\mathbf{y}$, we must evaluate the marginal distribution of $\mathcal{Y}$, which is the integral of $f_{\mathcal{Y},\mathcal{U}}(\mathbf{y},\mathbf{u})$ with respect to $\mathbf{u}$.

This is much simpler if we first rewrite the *penalized sum of squared residuals*, $\left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}
-\mathbf{Z}\boldsymbol{\Lambda}_\theta\mathbf{u}\right\|^2+
\left\|\mathbf{u}\right\|^2$, in @eq:YUjointdensity, which is a quadratic form in $\mathbf{u}$, to isolate the dependence on $\mathbf{u}$
$$
\begin{aligned}
  r^2_\theta(\mathbf{u},\boldsymbol{\beta})
  &=
  \|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}-\mathbf{Z}\boldsymbol{\Lambda}_\theta\mathbf{u}\|^2+\|\mathbf{u}\|^2 \\
  &=
  \left\|
    \begin{bmatrix}
      \mathbf{Z}\boldsymbol{\Lambda}_\theta & \mathbf{X} & \mathbf{y} \\
     -\mathbf{I}_q & \mathbf{0} & \mathbf{0}
    \end{bmatrix}
    \begin{bmatrix}
     -\mathbf{u} \\
     -\boldsymbol{\beta} \\
     1
    \end{bmatrix}
  \right\|^2 \\
  &=
    \begin{bmatrix}
     -\mathbf{u'} &
     -\boldsymbol{\beta}' &
      1
    \end{bmatrix}
    \begin{bmatrix}
      \boldsymbol{\Lambda}'\mathbf{Z}'\mathbf{Z}\boldsymbol{\Lambda}+\mathbf{I} & \boldsymbol{\Lambda}'\mathbf{Z}'\mathbf{X} & \boldsymbol{\Lambda}'\mathbf{Z}'\mathbf{y} \\
      \mathbf{X}'\mathbf{Z}\boldsymbol{\Lambda} & \mathbf{X}'\mathbf{X} & \mathbf{X}'\mathbf{y} \\
      \mathbf{y}'\mathbf{Z}\boldsymbol{\Lambda} & \mathbf{y}'\mathbf{X} & \mathbf{y}'\mathbf{y}
    \end{bmatrix}
    \begin{bmatrix}
     -\mathbf{u} \\
     -\boldsymbol{\beta} \\
      1
    \end{bmatrix} \\
     &=
    \begin{bmatrix}
     -\mathbf{u'} &
     -\boldsymbol{\beta'} &
      1
    \end{bmatrix}
    \begin{bmatrix}
      \mathbf{R}_{ZZ}' & \mathbf{0} & \mathbf{0} \\
      \mathbf{R}_{ZX}' & \mathbf{R}_{XX}' & \mathbf{0} \\
      \mathbf{r}_{Zy}' & \mathbf{r}_{Xy}' & r_{yy}
    \end{bmatrix}
    \begin{bmatrix}
      \mathbf{R}_{ZZ} & \mathbf{R}_{ZX} & \mathbf{r}_{Zy} \\
      \mathbf{0} & \mathbf{R}_{XX} & \mathbf{r}_{Xy} \\
      \mathbf{0} & \mathbf{0} & r_{yy}
    \end{bmatrix}
    \begin{bmatrix}
     -\mathbf{u} \\
     -\boldsymbol{\beta} \\
      1
    \end{bmatrix}\\
  &= \left\|
    \begin{bmatrix}
      \mathbf{R}_{ZZ} & \mathbf{R}_{ZX} & \mathbf{r}_{Zy}\\
      \mathbf{0} & \mathbf{R}_{XX}' & \mathbf{r}_{Xy}\\
      \mathbf{0} & \mathbf{0} & r_{yy}
    \end{bmatrix}
    \begin{bmatrix}
     -\mathbf{u} \\
     -\boldsymbol{\beta} \\
      1
    \end{bmatrix}
    \right\|^2\\
  &= \| \mathbf{r}_{Zy}-\mathbf{R}_{ZX}\boldsymbol{\beta}-\mathbf{R}_{ZZ}\mathbf{u} \|^2 +
     \| \mathbf{r}_{Xy}-\mathbf{R}_{XX}\boldsymbol{\beta}\|^2 + r_{yy}^2 ,
  \end{aligned}
$$ {#eq:penalized-rss}
using the Cholesky factor of the blocked matrix,
$$
\boldsymbol{\Omega}_\theta=
  \begin{bmatrix}
    \boldsymbol{\Lambda}_\theta'\mathbf{Z'Z}\boldsymbol{\Lambda}_\theta+\mathbf{I} & 
    \boldsymbol{\Lambda}_\theta'\mathbf{Z'X} & \boldsymbol{\Lambda}_\theta'\mathbf{Z'y} \\
    \mathbf{X'Z}\boldsymbol{\Lambda}_\theta & \mathbf{X'X} & \mathbf{X'y} \\
    \mathbf{y'Z}\boldsymbol{\Lambda}_\theta & \mathbf{y'X} & \mathbf{y'y}
  \end{bmatrix} =
  \begin{bmatrix}
    \mathbf{R}_{ZZ}' & \mathbf{0} & \mathbf{0} \\
    \mathbf{R}_{ZX}' & \mathbf{R}'_{XX} & \mathbf{0} \\
    \mathbf{r}_{Zy}' & \mathbf{r}'_{Xy} & r_{yy}
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{R}_{ZZ} & \mathbf{R}_{ZX} & \mathbf{r}_{Zy} \\
    \mathbf{0} & \mathbf{R}_{XX} & \mathbf{r}_{Xy} \\
    \mathbf{0} & \mathbf{0} & r_{yy}
  \end{bmatrix} .
$$ {#eq:bigCholfac}

Note that the block in the upper left, $\boldsymbol{\Lambda}_\theta'\mathbf{Z'Z}\boldsymbol{\Lambda}_\theta+\mathbf{I}$, is positive definite even when $\boldsymbol{\Lambda}_\theta$ is singular, because
$$
\mathbf{u}'\left(\boldsymbol{\Lambda}_\theta'\mathbf{Z'Z}\boldsymbol{\Lambda}_\theta+\mathbf{I}\right)\mathbf{u} = \left\|\mathbf{Z}\boldsymbol{\Lambda}_\theta\mathbf{u}\right\|^2
+\left\|\mathbf{u}\right\|^2
$$ {#eq:Cholfacupperleft}
and the first term is non-negative while the second is positive if $\mathbf{u}\ne\mathbf{0}$.

Thus $\mathbf{R}_{ZZ}$, with positive diagonal elements, can be evaluated and its determinant, $\left|\mathbf{R}_{ZZ}\right|$, is positive.
This determinant appears in the marginal density of $\mathcal{Y}$, from which the likelihood of the parameters is evaluated.

To evaluate the likelihood,
$$
L(\boldsymbol{\theta},\boldsymbol{\beta},\sigma|\mathbf{y}) = \int_\mathbf{u} f_{\mathcal{Y},\mathcal{U}}(\mathbf{y},\mathbf{u})\, d\mathbf{u}
$$ {#eq:likelihood-abstract}
we isolate the part of the joint density that depends on $\mathbf{u}$ and perform a change of variable
$$
\mathbf{v}=\mathbf{R}_{ZZ}\mathbf{u}+\mathbf{R}_{ZX}\boldsymbol{\beta}-\mathbf{r}_{Zy} .
$$ {#eq:u-system}
From the properties of the multivariate Gaussian distribution
$$
\begin{aligned}
  \int_{\mathbf{u}}\frac{1}{(2\pi\sigma^2)^{q/2}}
    \exp\left(-\frac{\|\mathbf{R}_{ZZ}\mathbf{u}+\mathbf{R}_{ZX}\boldsymbol{\beta}-\mathbf{r}_{Zy}\|^2}{2\sigma^2}\right)
    \,d\mathbf{u}
  &= \int_{\mathbf{v}}\frac{1}{(2\pi\sigma^2)^{q/2}}
    \exp\left(-\frac{\|\mathbf{v}\|^2}{2\sigma^2}\right)|\mathbf{R}_{ZZ}|^{-1}\,d\mathbf{v}\\
  &=|\mathbf{R}_{ZZ}|^{-1}
\end{aligned}
$$ {#eq:likelihood-integral}
from which we obtain the likelihood as
$$
  L(\boldsymbol{\theta},\boldsymbol{\beta},\sigma;\mathbf{y})=
  \frac{|\mathbf{R}_{ZZ}|^{-1}}{(2\pi\sigma^2)^{n/2}}
  \exp\left(-\frac{r_{yy}^2 + \|\mathbf{R}_{XX}(\boldsymbol{\beta}-\widehat{\boldsymbol{\beta}})\|^2}{2\sigma^2}\right) ,
$$ {#eq:likelihood}
where the conditional estimate, $\widehat{\boldsymbol{\beta}}$, given $\boldsymbol{\theta}$, satisfies
$$
\mathbf{R}_{XX}\widehat{\boldsymbol{\beta}} = \mathbf{r}_{Xy} .
$$

Setting $\boldsymbol{\beta}=\widehat{\boldsymbol{\beta}}$
and taking the logarithm provides the estimate of $\sigma^2$,
given $\boldsymbol{\theta}$, as
$$
\widehat{\sigma^2}=\frac{r_\mathbf{yy}^2}{n}
$$ {#eq:sigma-hat}
which gives the *profiled log-likelihood*,
$\ell(\boldsymbol{\theta}|\mathbf{y})=\log L(\boldsymbol{\theta},\widehat{\boldsymbol{\beta}},\widehat{\sigma})$,
on the deviance scale, as
$$
-2\ell(\boldsymbol{\theta}|\mathbf{y})=2\log(|\mathbf{R}_{ZZ}|) +
    n\left(1+\log\left(\frac{2\pi r_{yy}^2(\boldsymbol{\theta})}{n}\right)\right)
$$ {#eq:profiled-log-likelihood}

One of the interesting aspects of this formulation is that it is not necessary to solve for the conditional estimate of $\boldsymbol{\beta}$ or the conditional modes of the random effects when evaluating the log-likelihood.
The two values needed for the log-likelihood evaluation, $2\log(|\mathbf{R}_{ZZ}|)$ and $r_\mathbf{yy}^2$, are obtained directly from the diagonal elements of the Cholesky factor.

Furthermore, $\boldsymbol{\Omega}_{\theta}$ and, from that, the Cholesky factor, $\mathbf{R}_{\theta}$, and the objective to be optimized can be evaluated for a given value of $\boldsymbol{\theta}$ from
$$
\mathbf{A} = \begin{bmatrix}
\mathbf{Z}^\prime\mathbf{Z} & \mathbf{Z}^\prime\mathbf{X} & \mathbf{Z}^\prime\mathbf{y} \\
\mathbf{X}^\prime\mathbf{Z} & \mathbf{X}^\prime\mathbf{X} & \mathbf{X}^\prime\mathbf{y} \\
\mathbf{y}^\prime\mathbf{Z} & \mathbf{y}^\prime\mathbf{X} & \mathbf{y}^\prime\mathbf{y}
\end{bmatrix}
$$ {#eq:A}
and $\boldsymbol{\Lambda}_{\theta}$.

In the `MixedModels` package the `LinearMixedModel` struct contains a symmetric blocked array in the `A` field and a similarly structured lower-triangular blocked array in the `L` field.
Evaluation of the objective simply involves updating the template matrices, $\lambda_i, i=1,\dots,k$ in the `ReMat` structures then updating `L` from `A` and the $\lambda_i$.

## The REML criterion {#sec:REML}

The so-called REML estimates of variance components are often preferred to the maximum likelihood estimates.
("REML" can be considered to be an acronym for "restricted" or "residual" maximum likelihood, although neither term is completely accurate because these estimates do not maximize a likelihood.)
We can motivate the use of the REML criterion by considering a linear regression model, 
$$
  \mathcal{Y}\sim\mathcal{N}(\mathbf{X}\boldsymbol{\beta},\sigma^2\mathbf{I}_n),
$$ {#eq:20}
in which we typically estimate $\sigma^2$ as
$$
  \widehat{\sigma^2_R}=\frac{\|\mathbf{y}-\mathbf{X}\widehat{\boldsymbol{\beta}}\|^2}{n-p}
$$ {#eq:21}
even though the maximum likelihood estimate of $\sigma^2$ is
$$
  \widehat{\sigma^2_{L}}=\frac{\|\mathbf{y}-\vec
    X\widehat{\boldsymbol{\beta}}\|^2}{n} .
$$ {#eq:22}

The argument for preferring $\widehat{\sigma^2_R}$ to $\widehat{\sigma^2_{L}}$ as an estimate of $\sigma^2$ is that the numerator in both estimates is the sum of squared residuals at $\widehat{\boldsymbol{\beta}}$ and, although the residual vector, $\mathbf{y}-\mathbf{X}\widehat{\boldsymbol{\beta}}$, is an $n$-dimensional vector, it satisfies $p$ linearly independent constraints, $\mathbf{X}'(\mathbf{y}-\mathbf{X}\widehat{\boldsymbol{\beta}})=\mathbf{0}$.
That is, the residual at $\widehat{\boldsymbol{\beta}}$ is the projection of the observed response vector, $\mathbf{y}$, into an $(n-p)$-dimensional linear subspace of the $n$-dimensional response space.
The estimate $\widehat{\sigma^2_R}$ takes into account the fact that $\sigma^2$ is estimated from residuals that have only $n-p$ *degrees of freedom*.

Another argument often put forward for REML estimation is that $\widehat{\sigma^2_R}$ is an *unbiased* estimate of $\sigma^2$, in the sense that the expected value of the estimator is equal to the value of the parameter.
However, determining the expected value of an estimator involves integrating with respect to the density of the estimator and we have seen that densities of estimators of variances will be skewed, often highly skewed.
It is not clear why we should be interested in the expected value of a highly skewed estimator.
If we were to transform to a more symmetric scale, such as the estimator of the standard deviation or the estimator of the logarithm of the standard deviation, the REML estimator would no longer be unbiased.
Furthermore, this property of unbiasedness of variance estimators does not generalize from the linear regression model to linear mixed models.
This is all to say that the distinction between REML and ML estimates of variances and variance components is probably less important than many people believe.

Nevertheless it is worthwhile seeing how the computational techniques described in this chapter apply to the REML criterion because the REML parameter estimates $\widehat{\boldsymbol{\theta}}_R$ and $\widehat{\sigma_R^2}$ for a linear mixed model have the property that they would specialize to $\widehat{\sigma^2_R}$ from #eq:21 for a linear regression model, as seen in @sec:Dyestuff2LMM.

Although not usually derived in this way, the REML criterion (on the deviance scale) can be expressed as
$$
  d_R(\boldsymbol{\theta},\sigma|\mathbf{y})=-2\log
  \int_{\mathbb{R}^p}L(\boldsymbol{\theta},\boldsymbol{\beta},\sigma|\mathbf{y})\,d\boldsymbol{\beta} .
$$ {#eq:23}
The REML estimates $\widehat{\boldsymbol{\theta}}_R$ and $\widehat{\sigma_R^2}$
minimize $d_R(\boldsymbol{\theta},\sigma|\mathbf{y})$.

To evaluate this integral we form an expansion, similar to @eq:likelihood, of $r^2_{\theta,\beta}$ about $\widehat{\boldsymbol{\beta}}_\theta$
$$
  r^2_{\theta,\beta}=r^2_\theta+\|\mathbf{R}_{XX}(\boldsymbol{\beta}-\widehat{\boldsymbol{\beta}}_\theta)\|^2 .
$$ {#eq:rsqbetathetaexp}
from which we can derive
$$
  \int_{\mathbb{R}^p}\frac{\exp\left(-\frac{r^2_{\theta,\beta}}{2\sigma^2}\right)}
  {(2\pi\sigma^2)^{n/2}|\mathbf{R}_{ZZ}|} \,d\boldsymbol{\beta}=
  \frac{\exp\left(-\frac{r^2_\theta}{2\sigma^2}\right)}
  {(2\pi\sigma^2)^{(n-p)/2}|\mathbf{R}_{ZZ}||\mathbf{R}_X|}
$$ {#eq:betaintegral}
corresponding to a REML criterion on the deviance scale of
$$
  d_R(\boldsymbol{\theta},\sigma|\mathbf{y})=(n-p)\log(2\pi\sigma^2)+
  2\log\left(|\mathbf{R}_{ZZ}||\mathbf{R}_X|\right)+\frac{r^2_\theta}{\sigma^2} .
$$ {#eq:REMLdev}
Plugging in the conditional REML estimate, $\widehat{\sigma^2}_R=r^2_\theta/(n-p)$, provides the profiled REML criterion
$$
  \tilde{d}_R(\boldsymbol{\theta}|\mathbf{y})=
  2\log\left(|\mathbf{R}_{ZZ}||\mathbf{R}_X|\right)+(n-p)
  \left[1+\log\left(\frac{2\pi r^2_\theta}{n-p}\right)\right].
$$ {#eq:24}

The REML estimate of $\boldsymbol{\theta}$ is
$$
  \widehat{\boldsymbol{\theta}}_R=\arg\min_{\boldsymbol{\theta}}\tilde{d}_R(\boldsymbol{\theta}|\mathbf{y}) ,
$$ {#eq:31}
and the REML estimate of $\sigma^2$ is the conditional REML estimate of $\sigma^2$ at $\widehat{\boldsymbol{\theta}}_R$,
$$ \widehat{\sigma^2_R}=r^2_{\widehat\theta_R}/(n-p) . $$ {#eq:REMLsigmasq}
It is not entirely clear how one would define a "REML estimate" of $\boldsymbol{\beta}$
because the REML criterion, $d_R(\boldsymbol{\theta},\sigma|\mathbf{y})$, defined in @eq:REMLdev, does not depend on $\boldsymbol{\beta}$.
However, it is customary (and not unreasonable) to use
$\widehat{\boldsymbol{\beta}}_R=\widehat{\boldsymbol{\beta}}_{\widehat{\boldsymbol{\theta}}_R}$ as
the REML estimate of $\boldsymbol{\beta}$.