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

An alternative is to write the model in terms of the $n$-dimensional *response vector*, $\mathbf{y}$, an $n\times p$ *model matrix*, $\mathbf{X}, and a $p$-dimensional *coefficient vector*, $\boldsymbol{\beta}$, as
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
If $\sigma^2=0$ then all the probability is concentrated at a single point, $x=\mu$, and we no longer have a probability density, in the usual way of thinking of it.
Such a distribution is said to be [degenerate](https://en.wikipedia.org/wiki/Degenerate_distribution).

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
It also implies that there are "matrix square roots" of $\boldsymbol\Sigma$ in the sense that there are matrices $\mathbf{A}$ such that $\mathbf{A}\mathbf{A}'=\boldsymbol{\Sigma}$.
(The reason for writing $\mathbf{A}\mathbf{A}'$ and not simply the square of $\mathbf{A}$ is that $\mathbf{A}$ is not required to be symmetric but $\mathbf{A}\mathbf{A}'$ will be symmetric, even in $\mathbf{A}$ is not.)

One such "square root" of a positive definite $\boldsymbol\Sigma$ is the [Cholesky factor](https://en.wikipedia.org/wiki/Cholesky_decomposition), which is an $n\times n$ upper-triangular matrix, $\mathbf{R}$, such that
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
Requiring that the inverse of a matrix must be evaluated to solve a linear system is like saying that a quotient, $a/b$, must be evalated by calculating, $b^{-1}$, the reciprocal of $b$, then evaluating $b^{-1}a$, instead of evaluating the quotient directly.

In a derivation we may write an expression like $\mathbf{L}^{-1}\mathbf{b}$ but the evaluation is performed like @eq:forwardsolve.

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
Usually the log-likelihood is easier to optimize, either algebraically or numerically, than the likelihood itself.

The [maximum likelihood](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation) estimates of the parameters are, as the name implies, the values of $\boldsymbol{\beta}$ and $\sigma^2$ that maximize the expression on the right of @eq:linmodlikelihood .
Because the logarithm is a [monotone increasing](https://en.wikipedia.org/wiki/Monotonic_function) function, the maximum likelihood estimates will also maximize the *log-likelihood*

$$
\begin{align}
\ell(\boldsymbol{\beta},\sigma^2;\mathbf{y})
&=\log L(\boldsymbol{\beta},\sigma^2;\mathbf{y})\\
&=-\frac{n}{2}\log(2\pi\sigma^2)-\frac{\left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2}{2\sigma^2}
\end{align}
$$ {#eq:linmodloglike}

To avoid the negative signs and the factors of 2 in the denominator, we often convert the log-likelihood to the *deviance scale*, which is negative twice the log-likelihood
$$
\begin{align}
d(\boldsymbol{\beta},\sigma^2;\mathbf{y})
&=-2\ell(\boldsymbol{\beta},\sigma^2; \mathbf{y})\\
&=n\log(2\pi\sigma^2)+\frac{\left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2}{2\sigma^2}
\end{align}
$$ {#eq:linmoddevscale}
Because of the negative sign, the maximum likelihood estimates are those that *minimize* $d(\boldsymbol{\beta},\sigma^2;\mathbf{y})$.

(The term *deviance scale* is used for $d(\boldsymbol{\beta},\sigma^2;\mathbf{y})$ rather than [deviance](https://en.wikipedia.org/wiki/Deviance_(statistics)) because the deviance involves an additive shift, which is a correction for the saturated model - see the link.
It is obvious what the saturated model should be for the linear model but not for the linear mixed model so, to avoid confusion, we refer to the log-likelihood on the deviance scale as the *objective*.)
