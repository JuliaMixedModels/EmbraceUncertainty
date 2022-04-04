# Computational Methods for Mixed Models {#sec:computational}

In this chapter we describe some of the details of the computational methods for fitting linear mixed models, as implemented in the *MixedModels* package, and the theoretical development of these methods.
We also provide the basis for later generalizations to models for non-Gaussian responses and to models in which the relationship between the conditional mean, $\mathbf{\mu}$, and the linear predictor, $\mathbf{\gamma}=\mathbf{X}\mathbf{\beta}+\mathbf{Z}\mathbf{b}= \mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\mathbf{\beta}$, is a nonlinear relationship.

This material is directed at those readers who wish to follow the theory and methodology of linear mixed models and how both can be extended to other forms of mixed models.
Readers who are less interested in the "how" and the "why" of fitting mixed models than in the results themselves should not feel obligated to master these details.

We begin by reviewing the definition of linear mixed-effects models and some of the basics of the computational methods, as given in @bates.maechler.etal:2015 .

## Definitions and basic results {#sec:defnLMM}

As described in @bates.maechler.etal:2015 , a linear mixed-effects model is based on two vector-valued random variables: the $q$-dimensional vector of random effects, $\mathcal{B}$, and the $n$-dimensional response vector, $\mathcal{Y}$. @eq:LMMdist defines the unconditional distribution of $\mathcal{B}$ and the conditional distribution of $\mathcal{Y}$, given $\mathcal{B}=\mathbf{b}$, as multivariate Gaussian distributions of the form
$$
\begin{aligned}
  (\mathcal{Y}|\mathcal{B}=\mathbf{b})&\sim\mathcal{N}(\mathbf{X}\mathbf{\beta}+\mathbf{Z}\mathbf{b},\sigma^2\mathbf{I})\\
  \mathcal{B}&\sim\mathcal{N}(\mathbf{0},\Sigma_\theta) .
\end{aligned}
$$

The $q\times q$, symmetric, variance-covariance matrix, $\mathrm{Var}(\mathcal{B})=\Sigma_\theta$, depends on the *variance-component parameter vector*, $\mathbf{\theta}$, and is *positive semidefinite*, which means that
$$
  \mathbf{b}'\Sigma_\theta\mathbf{b}\ge0,\quad\forall\,\mathbf{b}\ne\mathbf{0} .
$$ {#eq:posSemiDef}
(The symbol $\forall$ denotes "for all".)
The fact that $\Sigma_\theta$ is positive semidefinite does not guarantee that $\Sigma_\theta^{-1}$ exists.
We would need a stronger property, $\mathbf{b}'\Sigma_\theta\mathbf{b}>0,\,\forall\,\mathbf{b}\ne\mathbf{0}$, called *positive definiteness*, to ensure that $\Sigma_\theta^{-1}$ exists.

Many computational formulas for linear mixed models are written in terms of this *precision matrix*, $\Sigma_\theta^{-1}$.
Such formulas will become unstable as $\Sigma_\theta$ approaches singularity.
And it can do so.
It is a fact that singular (i.e. non-invertible) $\Sigma_\theta$ can and do occur in practice, as we have seen in some of the examples in earlier chapters.
Moreover, during the course of the numerical optimization by which the parameter estimates are determined, it is frequently the case that the deviance or the REML criterion will need to be evaluated at values of $\mathbf{\theta}$ that produce a singular $\Sigma_\theta$.
Because of this we will take care to use computational methods that can be applied even when $\Sigma_\theta$ is singular and are stable as $\Sigma_\theta$ approaches singularity.

As defined in @eq:relcovfac, a relative covariance factor, $\Lambda_\theta$, is any matrix that satisfies
$$
\Sigma_\theta=\sigma^2\Lambda_\theta\Lambda_\theta' .
$$
According to this definition, $\Sigma$ depends on both $\sigma$ and $\theta$, and we should write it as $\Sigma_{\sigma,\theta}$.
However, we will blur that distinction and continue to write $\text{Var}(\mathcal{B})=\Sigma_\theta$.
Another technicality is that the *common scale parameter*, $\sigma$, could, in theory, be zero.
We will show that in practice the only way for its estimate, $\widehat{\sigma}$, to be zero is for the fitted values from the fixed-effects only, $\mathbf{X}\widehat{\mathbf{\beta}}$, to be exactly equal to the observed data.
This occurs only with data that have been (incorrectly) simulated without error.
In practice we can safely assume that $\sigma>0$.
However, $\Lambda_\theta$, like $\Sigma_\theta$, can be singular.

The computational methods in the *MixedModels* package are based on $\Lambda_\theta$ and do not require evaluation of $\Sigma_\theta$.
In fact, $\Sigma_\theta$ is explicitly evaluated only at the converged parameter estimates.

The spherical random effects, $\mathcal{U}\sim\mathcal{N}(\mathbf{0},\sigma^2\mathbf{I}_q)$, determine $\mathcal{B}$ as
$$
  \mathcal{B}=\Lambda_\theta\mathcal{U} .
$$ {#eq:sphericalRE}
Although it may seem more intuitive to write $\mathcal{U}$ as a linear transformation of $\mathcal{B}$, we cannot do that when $\Lambda_\theta$ is singular, which is why @eq:sphericalRE is in the form shown.

We can easily verify that @eq:sphericalRE provides the desired distribution for $\mathcal{B}$.
As a linear transformation of a multivariate Gaussian random variable, $\mathcal{B}$ will also be multivariate Gaussian.
Its mean and variance-covariance matrix are straightforward to evaluate,
$$
  \mathrm{E}[\mathcal{B}] = \mathrm{E}[\Lambda_\theta\mathcal{U}] = \Lambda_\theta\mathrm{E}[\mathcal{U}]=\Lambda_\theta\mathbf{0}=\mathbf{0}
$$ {#eq:EB}
and
$$
\mathrm{Var}(\mathcal{B})
  \begin{aligned}[t]
    &=\mathrm{E}\left[(\mathcal{B}-\mathrm{E}[\mathcal{B}])
      (\mathcal{B}-\mathrm{E}[\mathcal{B}])'\right]
    =\mathrm{E}\left[\mathcal{B}\mathcal{B}'\right]\\
    &=\mathrm{E}\left[\Lambda_\theta\,\mathcal{U}\mathcal{U}'\Lambda_\theta'\right]
    =\Lambda_\theta\,\mathrm{E}[\mathcal{U}\mathcal{U}']\Lambda_\theta'
    =\Lambda_\theta\,\mathrm{Var}(\mathcal{U})\Lambda_\theta'\\
    &=\Lambda_\theta\,\sigma^2\mathbf{I}_q\,\Lambda_\theta'
    =\sigma^2\Lambda_\theta\Lambda_\theta'
    =\Sigma_\theta
  \end{aligned}
$$
and have the desired form.

Just as we concentrate on how $\mathbf{\theta}$ determines $\Lambda_\theta$, not $\Sigma_\theta$, we will concentrate on properties of $\mathcal{U}$ rather than $\mathcal{B}$.
In particular, we now define the model according to the distributions
$$
  \begin{aligned}
  (\mathcal{Y}|\mathcal{U}=\mathbf{u})&\sim\mathcal{N}(\mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\beta,\sigma^2\mathbf{I}_n)\\
  \mathcal{U}&\sim\mathcal{N}(\mathbf{0},\sigma^2\mathbf{I}_q) .
  \end{aligned}
$$ {#eq:condYgivenU}

To allow for extensions to other types of mixed models we distinguish between the *linear predictor* 
$$
  \mathbf{\gamma} = \mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\beta
$$ {#eq:linearpred}
and the *conditional mean* of $\mathcal{Y}$, given $\mathcal{U}=\mathbf{u}$, which is
$$
  \mathbf{\mu} = \mathrm{E}\left[\mathcal{Y}|\mathcal{U}=\mathbf{u}\right] .
$$ {#eq:conditionalMean}
For a linear mixed model $\mathbf{\mu}=\mathbf{\gamma}$.
In other forms of mixed models the conditional mean, $\mathbf{\mu}$, can be a nonlinear function of the linear predictor, $\mathbf{\gamma}$.
For some models the dimension of $\mathbf{\gamma}$ is a multiple of $n$, the dimension of $\mathbf{\mu}$ and $\mathbf{y}$, but for a linear mixed model the dimension of $\mathbf{\gamma}$ must be $n$.
Hence, the model matrix $\mathbf{Z}$ must be $n\times q$ and $\mathbf{X}$ must be $n\times p$.

## The conditional distribution of U given Y {#sec:conddistUgivenY}

In this chapter it will help to be able to distinguish between the observed response vector and an arbitrary value of $\mathcal{Y}$.
For this chapter only we will write the observed data vector as $\mathbf{y}_{\text{obs}}$, with the understanding that $\mathbf{y}$ without the subscript will refer to an arbitrary value of the random variable $\mathcal{Y}$.

The likelihood of the parameters, $\mathbf{\theta}$, $\mathbf{\beta}$, and $\sigma$, given the observed data, $\mathbf{y}_{\text{obs}}$, is the probability density of $\mathcal{Y}$, evaluated at $\mathbf{y}_{\text{obs}}$.
Although the numerical values of the probability density and the likelihood are identical, the interpretations of these functions are different.
In the density we consider the parameters to be fixed and the value of $\mathbf{y}$ as varying.
In the likelihood we consider $\mathbf{y}$ to be fixed at $\mathbf{y}_{\text{obs}}$ and the parameters, $\mathbf{\theta}$, $\mathbf{\beta}$ and $\sigma$, as varying.

The natural approach for evaluating the likelihood is to determine the marginal distribution of $\mathcal{Y}$, which in this case amounts to determining the marginal density of $\mathcal{Y}$, and evaluate that density at $\mathbf{y}_{\text{obs}}$.
To follow this course we would first determine the joint density of $\mathcal{U}$ and $\mathcal{Y}$, written $f_{\mathcal{U},\mathcal{Y}}(\mathbf{u},\mathbf{y})$, then integrate this density with respect to $\mathbf{u}$ to create the marginal density, $f_{\mathcal{Y}}(\mathbf{y})$, and finally evaluate this marginal density at $\mathbf{y}_{\text{obs}}$.

To allow for later generalizations we will change the order of these steps slightly.
We evaluate the joint density function, $f_{\mathcal{U},\mathcal{Y}}(\mathbf{u},\mathbf{y})$, at $\mathbf{y}_{\text{obs}}$, producing the *unnormalized conditional density*, $h(\mathbf{u})$.
We say that $h$ is "unnormalized" because the conditional density is a multiple of $h$ 
$$
  f_{\mathcal{U}|\mathcal{Y}}(\mathbf{u}|\mathbf{y}_{\text{obs}})=\frac
  {h(\mathbf{u})}{\int_{\mathbb{R}^q}h(\mathbf{u})\,d\mathbf{u}}  .
$$ {#eq:conddenUgivenY}
In some theoretical developments the normalizing constant, which is the integral in the denominator of an expression like @eq:conddenUgivenY, is not of interest.
Here it is of interest because the normalizing constant is exactly the likelihood that we wish to evaluate,
$$
  L(\mathbf{\theta},\mathbf{\beta},\sigma|\mathbf{y}_{\text{obs}}) =
  \int_{\mathbb{R}^q}h(\mathbf{u})\,d\mathbf{u} .
$$ {#eq:LMMlikelihood}

For a linear mixed model, where all the distributions of interest are
multivariate Gaussian and the conditional mean, $\mathbf{\mu}$, is a linear
function of both $\mathbf{u}$ and $\mathbf{\beta}$, the distinction between
evaluating the joint density at $\mathbf{y}_{\text{obs}}$ to produce
$h(\mathbf{u})$ then integrating with respect to $\mathbf{u}$, as opposed to
first integrating the joint density then evaluating at $\mathbf{y}_{\text{obs}}$, is not terribly important.
For other mixed models this distinction can be important.
In particular, generalized linear mixed models, described in @sec:GLMMbinomial, are often used to model a discrete response, such as a binary response or a count, leading to a joint distribution for $\mathcal{Y}$ and $\mathcal{U}$ that is discrete with respect to one variable, $\mathbf{y}$, and continuous with respect to the other, $\mathbf{u}$.
In such cases there isn't a joint density for $\mathcal{Y}$ and $\mathcal{U}$. The necessary distribution theory for general $\mathbf{y}$ and $\mathbf{u}$ is well-defined but somewhat awkward to describe.
It is much easier to realize that we are only interested in the observed response vector, $\mathbf{y}_{\text{obs}}$, not some arbitrary value of $\mathbf{y}$, so we can concentrate on the conditional distribution of $\mathcal{U}$ given $\mathcal{Y}=\mathbf{y}_{\text{obs}}$. For all the mixed models we will consider, the conditional distribution, $(\mathcal{U}|\mathcal{Y}=\mathbf{y}_{\text{obs}})$, is continuous and both the conditional density, $f_{\mathcal{U}|\mathcal{Y}}(\mathbf{u}|\mathbf{y}_{\text{obs}})$, and its unnormalized form, $h(\mathbf{u})$, are well-defined.

## Integrating h in the linear mixed model {#sec:IntegratingH}

The integral defining the likelihood in @eq:LMMlikelihood has a closed form in the case of a linear mixed model but not for some of the more general forms of mixed models.
To motivate methods for approximating the likelihood in more general situations, we describe in some detail how the integral can be evaluated using a blocked Cholesky factor, $\mathbf{L}_\theta$, and the conditional mode,
$$ 
  \tilde{\mathbf{u}}=\arg\max_{\mathbf{u}} f_{\mathcal{U}|\mathcal{Y}}(\mathbf{u}|\mathbf{y}_{\text{obs}})=
  \arg\max_{\mathbf{u}} h(\mathbf{u}) = \arg\max_{\mathbf{u}}
  f_{\mathcal{Y}|\mathcal{U}}(\mathbf{y}_{\text{obs}}|\mathbf{u})\,f_{\mathcal{U}}(\mathbf{u}).
$$ {#eq:condMode}
The notation $\arg\max_{\mathbf{u}}$ means that $\tilde{\mathbf{u}}$ is the value of $\mathbf{u}$ that maximizes the expression that follows.

In general, the *mode* of a continuous distribution is the value of the random variable that maximizes the density.
The value $\tilde{\mathbf{u}}$ is called the conditional mode of $\mathbf{u}$, given $\mathcal{Y}=\mathbf{y}_{\text{obs}}$, because $\tilde{\mathbf{u}}$ maximizes the conditional density of $\mathcal{U}$ given $\mathcal{Y}=\mathbf{y}_{\text{obs}}$.
The location of the maximum can be determined by maximizing the unnormalized
conditional density because $h(\mathbf{u})$ is just a constant multiple of
$f_{\mathcal{U}|\mathcal{Y}}(\mathbf{u}|\mathbf{y}_{\text{obs}})$.
The last part of (eqn. @eq:condMode) is simply a re-expression of $h(\mathbf{u})$ as the product of $f_{\mathcal{Y}|\mathcal{U}}(\mathbf{y}_{\text{obs}}|\mathbf{u})$ and $f_{\mathcal{U}}(\mathbf{u})$.
For a linear mixed model these densities are
$$
\begin{aligned}
  \label{eq:densYgivenUandU}
  f_{\mathcal{Y}|\mathcal{U}}(\mathbf{y}|\mathbf{u})&=
  \frac{1}{\left(2\pi\sigma^2\right)^{n/2}}
  \exp\left(-\frac{\left\|\mathbf{y}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2}{2\sigma^2}\right)\\
  f_{\mathcal{U}}(\mathbf{u})&=
  \frac{1}{\left(2\pi\sigma^2\right)^{q/2}}\exp\left(-\frac{\|\mathbf{u}\|^2}
    {2\sigma^2}\right)
\end{aligned}
$$
with product
$$
  h(\mathbf{u})=\frac{1}{\left(2\pi\sigma^2\right)^{(n+q)/2}}
  \exp\left(-\frac{\left\|\mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+\|\mathbf{u}\|^2}{2\sigma^2}\right) .
$$ {#eq:hudef}
On the deviance scale we have
$$
  -2\log\left(h(\mathbf{u})\right)=(n+q)\log(2\pi\sigma^2)
  +\frac{\left\|\mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+\|\mathbf{u}\|^2}{\sigma^2} .
$$ {#eq:devh}
Because @eq:devh describes the negative log density, $\tilde{\mathbf{u}}$ will be the value of $\mathbf{u}$ that minimizes the expression on the right hand side of @eq:devh.

The only part of the right hand side of @eq:devh that depends on $\mathbf{u}$ is the numerator of the second term.
Thus 
$$
  \tilde{\mathbf{u}}=\arg\min_{\mathbf{u}} \left\|
    \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+
  \|\mathbf{u}\|^2.
$$ {#eq:PLSsol}
The expression to be minimized, called the *objective function*, is described as a *penalized residual sum of squares* (PRSS) and the minimizer, $\tilde{\mathbf{u}}$, is called the *penalized least squares* (PLS) solution.
They are given these names because the first term in the objective, $\left\| \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2$, is a sum of squared residuals, and the second term, $\|\mathbf{u}\|^2$, is a penalty on the length, $\|\mathbf{u}\|$, of $\mathbf{u}$.
Larger values of $\mathbf{u}$ (in the sense of greater lengths as vectors) incur a higher penalty.

The PRSS criterion determining the conditional mode balances fidelity to the observed data (i.e. producing a small residual sum of squares) against simplicity of the model (small $\|\mathbf{u}\|$).
We refer to this type of criterion as a smoothing objective, in the sense that it seeks to smooth out the fitted response by reducing model complexity while still retaining reasonable fidelity to the observed data.

For the purpose of evaluating the likelihood we will regard the PRSS criterion as a function of the parameters, $\theta$ and $\beta$, given the data, $\mathbf{y}_{\text{obs}}$, and write its minimum value as
$$
  r^2_{\theta,\beta}=\min_{\mathbf{u}} \left\|
    \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+ \|\mathbf{u}\|^2.
$$ {#eq:r2thetabeta}
Notice that $\mathbf{\beta}$ enters the right hand side of @eq:r2thetabeta only through the linear predictor expression.
We will see that $\tilde{\mathbf{u}}$ can be determined by a direct (i.e. non-iterative) calculation and, in fact, we can minimize the PRSS criterion with respect to $\mathbf{u}$ and $\mathbf{\beta}$ simultaneously without iterating.
We write this minimum value as
$$
  r^2_\theta=\min_{\mathbf{u},\mathbf{\beta}} \left\|
    \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+ \|\mathbf{u}\|^2.
$$ {#eq:r2theta}
The value of $\mathbf{\beta}$ at the minimum is called the conditional estimate of $\mathbf{\beta}$ given $\mathbf{\theta}$, written $\widehat{\mathbf{\beta}}_\theta$.

## Determining the PLS solutions

As described in @bates.maechler.etal:2015 the profiled log-likelihood or profiled REML criterion for an LMM can be evaluated from the minimum of a particular penalized least squares (PLS) problem and determinants of Cholesky factors that are used to solve that problem.

This is based on viewing the relative covariance factor $\mathbf{\Lambda}_\theta$ as generating the random effects vector, $\mathcal{B}$, with distribution $\mathcal{N}(\mathbf{0},\Sigma)$, from a *spherical* random effects vector, $\mathcal{U}$, as
$$
\mathcal{B} = \mathbf{\Lambda}_{\mathbf{\theta}} \mathcal{U}\quad\mathrm{where}\quad\mathcal{U}\sim\mathcal{N}(\mathbf{0},\sigma^2\mathbf{I}_q) .
$$ {#eq:spherical-def}

The joint density of $\mathcal{Y}$ and $\mathcal{U}$, which is the product of the conditional density,
$f_{\mathcal{Y}|\mathcal{U}=\mathbf{u}}(\mathbf{y}|\mathbf{u})$, and
the unconditional density, $f_{\mathcal{U}}(\mathbf{u})$, becomes
$$
f_{\mathcal{Y},\mathcal{U}}(\mathbf{y,u})=\frac{1}{(2\pi\sigma^2)^{(n+q)/2}}\exp\left(\frac{-r^2_\mathbf{\theta}(\mathbf{u},\mathbf{\beta})}{2\sigma^2}\right)
$$ {#eq:joint-density}
where the penalized sum of squared residuals, $r^2_\mathbf{\theta}(\mathbf{u},\mathbf{\beta})$, is
$$
\begin{aligned}
  r^2_\mathbf{\theta}(\mathbf{u},\mathbf{\beta})
  &=
  \|\mathbf{y}-\mathbf{X\beta}-\mathbf{Z\Lambda_\theta u}\|^2+\|\mathbf{u}\|^2 \\
  &=
  \left\|
    \begin{bmatrix}
      \mathbf{Z\Lambda} & \mathbf{X} & \mathbf{y} \\
     -\mathbf{I}_q & \mathbf{0} & \mathbf{0}
    \end{bmatrix}
    \begin{bmatrix}
     -\mathbf{u} \\
     -\mathbf{\beta} \\
     1
    \end{bmatrix}
  \right\|^2 \\
  &=
    \begin{bmatrix}
     -\mathbf{u^\prime} &
     -\mathbf{\beta^\prime} &
      1
    \end{bmatrix}
    \begin{bmatrix}
      \mathbf{\Lambda}^\prime\mathbf{Z}^\prime\mathbf{Z\Lambda}+\mathbf{I} & \mathbf{\Lambda}^\prime\mathbf{Z}^\prime\mathbf{X} & \mathbf{\Lambda}^\prime\mathbf{Z}^\prime\mathbf{y} \\
      \mathbf{X}^\prime\mathbf{Z\Lambda} & \mathbf{X}^\prime\mathbf{X} & \mathbf{X}^\prime\mathbf{y} \\
      \mathbf{y}^\prime\mathbf{Z\Lambda} & \mathbf{y}^\prime\mathbf{X} & \mathbf{y}^\prime\mathbf{y}
    \end{bmatrix}
    \begin{bmatrix}
     -\mathbf{u} \\
     -\mathbf{\beta} \\
      1
    \end{bmatrix} \\
     &=
    \begin{bmatrix}
     -\mathbf{u^\prime} &
     -\mathbf{\beta^\prime} &
      1
    \end{bmatrix}
    \begin{bmatrix}
      \mathbf{L_{ZZ}} & \mathbf{0} & \mathbf{0} \\
      \mathbf{L_{XZ}} & \mathbf{L_{XX}} & \mathbf{0} \\
      \mathbf{l_{yZ}} & \mathbf{l_{yX}} & l_\mathbf{yy}
    \end{bmatrix}
    \begin{bmatrix}
      \mathbf{L_{ZZ}^\prime} & \mathbf{L_{XZ}^\prime} & \mathbf{l_{yZ}^\prime} \\
      \mathbf{0} & \mathbf{L_{XX}^\prime} & \mathbf{l_{yX}^\prime} \\
      \mathbf{0} & \mathbf{0} & l_\mathbf{yy}
    \end{bmatrix}
    \begin{bmatrix}
     -\mathbf{u} \\
     -\mathbf{\beta} \\
      1
    \end{bmatrix}\\
  &= \left\|
    \begin{bmatrix}
      \mathbf{L_{ZZ}^\prime} & \mathbf{L_{XZ}^\prime} & \mathbf{l_{yZ}^\prime} \\
      \mathbf{0} & \mathbf{L_{XX}^\prime} & \mathbf{l_{yX}^\prime} \\
      \mathbf{0} & \mathbf{0} & l_\mathbf{yy}
    \end{bmatrix}
    \begin{bmatrix}
     -\mathbf{u} \\
     -\mathbf{\beta} \\
      1
    \end{bmatrix}
    \right\|^2\\
  &= \| \mathbf{l_{yZ}^\prime}-\mathbf{L_{XZ}^\prime\beta}-\mathbf{L_{ZZ}\mathbf{u}} \|^2 +
     \| \mathbf{l_{yX}^\prime}-\mathbf{L_{XX}^\prime\beta}\|^2 + l_\mathbf{yy}^2 .
  \end{aligned}
$$ {#eq:penalized-rss}
using the blocked lower-triangular, sparse Cholesky factor,
$$
  \mathbf{\Omega_\theta}=
  \begin{bmatrix}
    \mathbf{\Lambda_\theta^\prime Z^\prime Z\Lambda_\theta+I} & \mathbf{\Lambda_\theta^\prime Z^\prime X} & \mathbf{\Lambda_\theta^\prime Z^\prime y} \\
    \mathbf{X^\prime Z\Lambda_\theta} & \mathbf{X^\prime X} & \mathbf{X^\prime y} \\
    \mathbf{y^\prime Z\Lambda_\theta} & \mathbf{y^\prime X} & \mathbf{y^\prime y}
  \end{bmatrix} =
  \begin{bmatrix}
    \mathbf{L_{ZZ}} & \mathbf{0} & \mathbf{0} \\
    \mathbf{L_{XZ}} & \mathbf{L_{XX}} & \mathbf{0} \\
    \mathbf{l_{yZ}} & \mathbf{l_{yX}} & l_\mathbf{yy}
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{L_{ZZ}^\prime} & \mathbf{L_{XZ}^\prime} & \mathbf{l_{yZ}^\prime} \\
    \mathbf{0} & \mathbf{L_{XX}^\prime} & \mathbf{l_{yX}^\prime} \\
    \mathbf{0} & \mathbf{0} & l_\mathbf{yy}
  \end{bmatrix} .
$$ {#eq:Omega}

From @eq:penalized-rss it is obvious that, for a fixed value of $\mathbf{\theta}$, the minimum
$r^2_\mathbf{\theta}(\mathbf{u},\mathbf{\beta})$ is $l_\mathbf{yy}^2$ and
the conditional estimate of $\mathbf{\beta}$ satisfies
$$
\mathbf{L_{XX}^\prime}\widehat{\mathbf{\beta}}(\mathbf{\theta})=\mathbf{l_{yX}^\prime} .
$$
The conditional mode, $\tilde{\mathbf{u}}$, of
$\mathcal{U}$ given $\mathcal{Y}=\mathbf{y}$, is the solution to
$$
\mathbf{L_{ZZ}^\prime}\tilde{\mathbf{u}}=\mathbf{l_{yZ}^\prime}-\mathbf{L_{XZ}^\prime}\widehat{\mathbf{\beta}} .
$$ {#eq:condmode}
(Technically, $\mathbf{\beta}$ and $\mathbf{\theta}$ in @eq:condmode are assumed known because this expression is a statement about distributions. 
In practice, the estimates, $\widehat{\mathbf{\theta}}$ and $\widehat{\beta}$, are plugged in when evaluating the conditional modes or "best linear unbiased predictors (BLUPs)" as they are sometimes called.)

The determinant, $|\mathbf{L_{ZZ}}|$, is the product of its diagonal elements, which must be positive.
Assuming that the fixed-effects model matrix, $\mathbf{X}$, has full column rank (this is checked and, if necessary, adjusted in a pre-processing step), $\mathbf{L_{XX}}$ also has positive diagonal elements.

To evaluate the likelihood,
$$
L(\mathbf{\theta},\mathbf{\beta},\sigma|\mathbf{y}) = \int_\mathbf{u} f_{\mathcal{Y},\mathcal{U}}(\mathbf{y},\mathbf{u})\, d\mathbf{u}
$$ {#eq:likelihood-abstract}
we isolate the part of the joint density that depends on $\mathbf{u}$ and perform a change of variable
$$
\mathbf{v}=\mathbf{L_{ZZ}^\prime u}+\mathbf{L_{XZ}^\prime\beta}-\mathbf{l_{yZ}^\prime} .
$$ {#eq:u-system}
From the properties of the multivariate Gaussian distribution
$$
\begin{aligned}
  \int_{\mathbf{u}}\frac{1}{(2\pi\sigma^2)^{q/2}}
    \exp\left(-\frac{\|\mathbf{L_{ZZ}^\prime u}+\mathbf{L_{XZ}^\prime\beta}-\mathbf{l_{yZ}^\prime}\|^2}{2\sigma^2}\right)
    \,d\mathbf{u}
  &= \int_{\mathbf{v}}\frac{1}{(2\pi\sigma^2)^{q/2}}
    \exp\left(-\frac{\|\mathbf{v}\|^2}{2\sigma^2}\right)|\mathbf{L_{ZZ}^\prime}|^{-1}\,d\mathbf{v}\\
  &=|\mathbf{L_{ZZ}}|^{-1}
\end{aligned}
$$ {#eq:likelihood-integral}
from which we obtain the likelihood as
$$
  L(\mathbf{\theta},\mathbf{\beta},\sigma)=
  \frac{|\mathbf{L_{ZZ}}|^{-1}}{(2\pi\sigma^2)^{n/2}}
  \exp\left(-\frac{l_\mathbf{yy}^2 + \|\mathbf{L_{XX}^\prime}(\mathbf{\beta}-\widehat{\mathbf{\beta}})\|^2}{2\sigma^2}\right) .
$$ {#eq:likelihood}
Setting $\mathbf{\beta}=\widehat{\mathbf{\beta}}$
and taking the logarithm provides the estimate of $\sigma^2$,
given $\mathbf{\theta}$, as
$$
\widehat{\sigma^2}=\frac{l_\mathbf{yy}^2}{n}
$$ {#eq:sigma-hat}
which gives the *profiled log-likelihood*,
$\ell(\mathbf{\theta}|\mathbf{y})=\log L(\mathbf{\theta},\widehat{\mathbf{\beta}},\widehat{\sigma})$,

as
$$
-2\ell(\mathbf{\theta}|\mathbf{y})=2\log(|\mathbf{L_{ZZ}}|) +
    n\left(1+\log\left(\frac{2\pi l_\mathbf{yy}^2(\mathbf{\theta})}{n}\right)\right)
$$ {#eq:profiled-log-likelihood}

This may seem complicated but, relative to other formulations of the model, it is remarkably simple.

One of the interesting aspects of this formulation is that it is not necessary to solve for the conditional estimate of $\mathbf{\beta}$ or the conditional modes of the random effects when evaluating the log-likelihood.
The two values needed for the log-likelihood evaluation, $2\log(|\mathbf{L}_{ZZ}|)$ and $l_\mathbf{yy}^2$, are obtained directly from the diagonal elements of the Cholesky factor.

Furthermore, $\mathbf{\Omega_\theta}$ and, from that, the Cholesky factor, $\mathbf{L_\theta}$, and the objective to be optimized can be evaluated for a given value of $\mathbf{\theta}$ from
$$
\mathbf{A} = \begin{bmatrix}
\mathbf{Z}^\prime\mathbf{Z} & \mathbf{Z}^\prime\mathbf{X} & \mathbf{Z}^\prime\mathbf{y} \\
\mathbf{X}^\prime\mathbf{Z} & \mathbf{X}^\prime\mathbf{X} & \mathbf{X}^\prime\mathbf{y} \\
\mathbf{y}^\prime\mathbf{Z} & \mathbf{y}^\prime\mathbf{X} & \mathbf{y}^\prime\mathbf{y}
\end{bmatrix}
$$ {#eq:A}
and $\mathbf{\Lambda_\theta}$.

In the `MixedModels` package the `LinearMixedModel` struct contains a symmetric blocked array in the `A` field and a similarly structured lower-triangular blocked array in the `L` field.
Evaluation of the objective simply involves updating the template matrices, $\lambda_i, i=1,\dots,k$ in the `ReMat` structures then updating `L` from `A` and the $\lambda_i$.

## The REML criterion {#sec:REML}

The so-called REML estimates of variance components are often preferred to the maximum likelihood estimates.
("REML" can be considered to be an acronym for "restricted" or "residual" maximum likelihood, although neither term is completely accurate because these estimates do not maximize a likelihood.)
We can motivate the use of the REML criterion by considering a linear regression model, 
$$
  \mathcal{Y}\sim\mathcal{N}(\mathbf{X}\mathbf{\beta},\sigma^2\mathbf{I}_n),
$$ {#eq:20}
in which we typically estimate $\sigma^2$ as
$$
  \widehat{\sigma^2_R}=\frac{\|\mathbf{y}_{\text{obs}}-\mathbf{X}\widehat{\mathbf{\beta}}\|^2}{n-p}
$$ {#eq:21}
even though the maximum likelihood estimate of $\sigma^2$ is
$$
  \widehat{\sigma^2_{L}}=\frac{\|\mathbf{y}_{\text{obs}}-\vec
    X\widehat{\mathbf{\beta}}\|^2}{n} .
$$ {#eq:22}

The argument for preferring $\widehat{\sigma^2_R}$ to $\widehat{\sigma^2_{L}}$ as an estimate of $\sigma^2$ is that the numerator in both estimates is the sum of squared residuals at $\widehat{\mathbf{\beta}}$ and, although the residual vector, $\mathbf{y}_{\text{obs}}-\mathbf{X}\widehat{\mathbf{\beta}}$, is an $n$-dimensional vector, it satisfies $p$ linearly independent constraints, $\mathbf{X}'(\mathbf{y}_{\text{obs}}-\mathbf{X}\widehat{\mathbf{\beta}})=\mathbf{0}$.
That is, the residual at $\widehat{\mathbf{\beta}}$ is the projection of the observed response vector, $\mathbf{y}_{\text{obs}}$, into an $(n-p)$-dimensional linear subspace of the $n$-dimensional response space.
The estimate $\widehat{\sigma^2_R}$ takes into account the fact that $\sigma^2$ is estimated from residuals that have only $n-p$ *degrees of freedom*.

Another argument often put forward for REML estimation is that $\widehat{\sigma^2_R}$ is an *unbiased* estimate of $\sigma^2$, in the sense that the expected value of the estimator is equal to the value of the parameter.
However, determining the expected value of an estimator involves integrating with respect to the density of the estimator and we have seen that densities of estimators of variances will be skewed, often highly skewed.
It is not clear why we should be interested in the expected value of a highly skewed estimator.
If we were to transform to a more symmetric scale, such as the estimator of the standard deviation or the estimator of the logarithm of the standard deviation, the REML estimator would no longer be unbiased.
Furthermore, this property of unbiasedness of variance estimators does not generalize from the linear regression model to linear mixed models.
This is all to say that the distinction between REML and ML estimates of variances and variance components is probably less important than many people believe.

Nevertheless it is worthwhile seeing how the computational techniques described in this chapter apply to the REML criterion because the REML parameter estimates $\widehat{\mathbf{\theta}}_R$ and $\widehat{\sigma_R^2}$ for a linear mixed model have the property that they would specialize to $\widehat{\sigma^2_R}$ from #eq:21 for a linear regression model, as seen in @sec:Dyestuff2LMM.

Although not usually derived in this way, the REML criterion (on the deviance scale) can be expressed as
$$
  d_R(\mathbf{\theta},\sigma|\mathbf{y}_{\text{obs}})=-2\log
  \int_{\mathbb{R}^p}L(\mathbf{\theta},\mathbf{\beta},\sigma|\mathbf{y}_{\text{obs}})\,d\mathbf{\beta} .
$$ {#eq:23}
The REML estimates $\widehat{\mathbf{\theta}}_R$ and $\widehat{\sigma_R^2}$
minimize $d_R(\mathbf{\theta},\sigma|\mathbf{y}_{\text{obs}})$.

To evaluate this integral we form an expansion, similar to @eq:PRSSwithL, of $r^2_{\theta,\beta}$ about $\widehat{\mathbf{\beta}}_\theta$
$$
  r^2_{\theta,\beta}=r^2_\theta+\|\mathbf{R}_X(\mathbf{\beta}-\widehat{\mathbf{\beta}}_\theta)\|^2 .
$$ {#eq:rsqbetathetaexp}
In the same way that @eq:PRSSwithL was used to simplify the integral in @eq:hintegral, we can derive
$$
  \int_{\mathbb{R}^p}\frac{\exp\left(-\frac{r^2_{\theta,\beta}}{2\sigma^2}\right)}
  {(2\pi\sigma^2)^{n/2}|\mathbf{L}_\theta|} \,d\mathbf{\beta}=
  \frac{\exp\left(-\frac{r^2_\theta}{2\sigma^2}\right)}
  {(2\pi\sigma^2)^{(n-p)/2}|\mathbf{L}_\theta||\mathbf{R}_X|}
$$ {#eq:betaintegral}
corresponding to a REML criterion on the deviance scale of
$$
  d_R(\mathbf{\theta},\sigma|\mathbf{y}_{\text{obs}})=(n-p)\log(2\pi\sigma^2)+
  2\log\left(|\mathbf{L}_\theta||\mathbf{R}_X|\right)+\frac{r^2_\theta}{\sigma^2} .
$$ {#eq:REMLdev}
Plugging in the conditional REML estimate, $\widehat{\sigma^2}_R=r^2_\theta/(n-p)$, provides the profiled REML criterion
$$
  \tilde{d}_R(\mathbf{\theta}|\mathbf{y}_{\text{obs}})=
  2\log\left(|\mathbf{L}_\theta||\mathbf{R}_X|\right)+(n-p)
  \left[1+\log\left(\frac{2\pi r^2_\theta}{n-p}\right)\right].
$$ {#eq:24}

The REML estimate of $\mathbf{\theta}$ is
$$
  \widehat{\mathbf{\theta}}_R=\arg\min_{\mathbf{\theta}}\tilde{d}_R(\mathbf{\theta}|\mathbf{y}_{\text{obs}}) ,
$$ {#eq:31}
and the REML estimate of $\sigma^2$ is the conditional REML estimate of $\sigma^2$ at $\widehat{\mathbf{\theta}}_R$,
$$ \widehat{\sigma^2_R}=r^2_{\widehat\theta_R}/(n-p) . $$ {#eq:REMLsigmasq}
It is not entirely clear how one would define a "REML estimate" of $\mathbf{\beta}$
because the REML criterion, $d_R(\mathbf{\theta},\sigma|\mathbf{y})$, defined in @eq:REMLdev, does not depend on $\mathbf{\beta}$.
However, it is customary (and not unreasonable) to use
$\widehat{\mathbf{\beta}}_R=\widehat{\mathbf{\beta}}_{\widehat{\mathbf{\theta}}_R}$ as
the REML estimate of $\mathbf{\beta}$.

## Step-by-step evaluation of the profiled deviance {#sec:stepByStep}

## Generalizing to other forms of mixed models {#sec:general}

In later chapters we cover the theory and practice of generalized linear mixed models (GLMMs), nonlinear mixed models (NLMMs) and generalized nonlinear mixed models (GNLMMs).
Because quite a bit of the theoretical and computational methodology covered in this chapter extends to those models we will cover the common aspects here.

### Descriptions of the model forms {#sec:modelForms}

We apply the name "generalized" to models in which the conditional distribution, $(\mathcal{Y}|\mathcal{U}=\mathbf{u})$, is not required to be Gaussian but does preserve some of the properties of the spherical Gaussian conditional distribution
$$
(\mathcal{Y}|\mathcal{U}=\mathbf{u})\sim\mathcal{N}
  (\mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\mathbf{\beta},\sigma^2\mathbf{I}_n)
$$
from the linear mixed model.
In particular, the components of $\mathcal{Y}$ are *conditionally independent*, given $\mathcal{U}=\mathbf{u}$.
Furthermore, $\mathbf{u}$ affects the distribution only through the conditional mean, which we will continue to write as $\mathbf{\mu}$, and it affects the conditional mean only through the linear predictor, $\mathbf{\gamma}=\mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\mathbf{\beta}$.

Typically we do not have $\mathbf{\mu}=\mathbf{\gamma}$, however.
The elements of the linear predictor, $\mathbf{\gamma}$, can be positive or negative or zero.
Theoretically they can take on any value between $-\infty$ and $\infty$.
But many distributional forms used in GLMMs put constraints on the value of the mean.
For example, the mean of a Bernoulli random variable, modeling a binary response, must be in the range $0<\mu<1$ and the mean of a Poisson random variable, modeling a count, must be positive.
To achieve these constraints we write the conditional mean, $\mathbf{\mu}$, as a transformation of the unbounded predictor, written $\mathbf{\eta}$.
For historical, and some theoretical, reasons the inverse of this transformation is called the *link function*, written
$$
  \mathbf{\eta}=\mathbf{g}(\mathbf{\mu}) ,
$$ {#eq:linkfunction}
and the transformation we want is called the *inverse link*, written $\mathbf{g}^{-1}$.

Both $\mathbf{g}$ and $\mathbf{g}^{-1}$ are determined by scalar functions, $g$ and $g^{-1}$, respectively, applied to the individual components of the vector argument.
That is, $\mathbf{\eta}$ must be $n$-dimensional and the vector-valued function $\mathbf{\mu}=\mathbf{g}^{-1}(\mathbf{\eta})$ is defined by the component functions $\mu_i=g^{-1}(\eta_i),\,i=1,\dots,n$.
Among other things, this means that the Jacobian matrix of the inverse link, $\frac{d\mathbf{\mu}}{d\mathbf{\eta}}$, will be diagonal.

Because the link function, $\mathbf{g}$, and the inverse link, $\mathbf{g}^{-1}$, are nonlinear functions (there would be no purpose in using a linear link function) many people use the terms "generalized linear mixed model" and "nonlinear mixed model" interchangeably.
We reserve the term "nonlinear mixed model" for the type of models used, for example, in pharmacokinetics and pharmacodynamics, where the conditional distribution is a spherical multivariate Gaussian
$$
  (\mathcal{Y}|\mathcal{U}=\mathbf{u})\sim\mathcal{N}(\mathbf{\mu}, \sigma^2\mathbf{I}_n)
$$ {#eq:NLMMconddist}
but
$\mathbf{\mu}$ depends nonlinearly on $\mathbf{\gamma}$.
For NLMMs the length of the linear predictor, $\mathbf{\gamma}$, is a multiple, $ns$, of $n$, the length of $\mathbf{\mu}$.

Like the map from $\mathbf{\eta}$ to $\mathbf{\mu}$, the map from $\mathbf{\gamma}$ to $\mathbf{\mu}$ has a "diagonal" property, which we now describe.
If we use $\mathbf{\gamma}$ to fill the columns of an $n\times s$ matrix, $\Gamma$, then $\mu_i$ depends only on the $i$th row of $\Gamma$.
In fact, $\mu_i$ is determined by a nonlinear model function, $f$, applied to the $i$ row of $\Gamma$.
Writing $\mathbf{\mu}=\mathbf{f}(\mathbf{\gamma})$ based on the component function $f$, we see that the Jacobian of $\mathbf{f}$, $\frac{d\mathbf{\mu}}{d\mathbf{\gamma}}$, will be the vertical concatenation of $s$ diagonal $n\times n$ matrices.

Because we will allow for generalized nonlinear mixed models (GNLMMs), in which the mapping from $\mathbf{\gamma}$ to $\mathbf{\mu}$ has the form
$$
  \mathbf{\gamma}\;\rightarrow\;\mathbf{\eta}\;\rightarrow\;\mathbf{\mu} ,
$$ {#eq:GammaEtaMu}
we will use @eq:GammaEtaMu in our definitions.

### Determining the conditional mode {#sec:conditionalmode}

For all these types of mixed models, the conditional distribution, $(\mathcal{U}|\mathcal{Y}=\mathbf{y}_{\text{obs}})$, is a continuous distribution for which we can determine the unscaled conditional density, $h(\mathbf{u})$.
As for linear mixed models, we define the conditional mode, $\tilde{\mathbf{u}}$, as the value that maximizes the unscaled conditional density.

Determining the conditional mode, $\tilde{\mathbf{u}}$, in a nonlinear mixed model is a penalized nonlinear least squares (PNLS) problem 
$$
  \tilde{\mathbf{u}}=\arg\min_{\mathbf{u}}\|\mathbf{y}_{\text{obs}}-\mathbf{\mu}\|^2+\|\mathbf{u}\|^2
$$ {#eq:PNLSprob}
which we solve by adapting the iterative techniques, such as the Gauss-Newton method [@bateswatts88:_nonlin Sect. 2.2.1], used for nonlinear least squares.
Starting at an initial value, $\mathbf{u}^{(0)}$, (the bracketed superscript denotes the iteration number) with conditional mean, $\mathbf{\mu}^{(0)}$, we determine an increment 
$\mathbf{\delta}^{(1)}$ by solving the penalized linear least squares
problem,
$$
  \mathbf{\delta}^{(1)}=\arg\min_{\mathbf{\delta}}\left\|
    \begin{bmatrix}
      \mathbf{y}_{\text{obs}}-\mathbf{\mu}^{(0)}\\
      \mathbf{0}-\mathbf{u}^{(0)}
    \end{bmatrix} -
    \begin{bmatrix}
      \mathbf{u}^{(0)}\\
      \mathbf{I}_q
    \end{bmatrix}\mathbf{\delta}\right\|^2$$ where $$\label{eq:Udef}
  \mathbf{u}^{(0)}=\left.\frac{d\mathbf{\mu}}{d\mathbf{u}}\right|_{\mathbf{u}^{(0)}} .
$$ {#eq:PNLSiter}
Naturally, we use the blocked Cholesky decomposition, $\mathbf{L}_\theta^{(0)}$, satisfying $$
  \mathbf{L}_\theta^{(0)}\left(\mathbf{L}_\theta^{(0)}\right)=
  \left[\left(\mathbf{u}^{(0)}\right)'\mathbf{u}^{(0)}+\mathbf{I}_q\right]
$$ {#eq:blockedCholPNLS}
to determine this increment.
The next iteration begins at
$$
  \mathbf{u}^{(1)} = \mathbf{u}^{(0)} + k \mathbf{\delta}^{(1)}
$$ {#eq:u1}
where $k$ is the step factor chosen, perhaps by step-halving [@bateswatts88:_nonlin Sect. 2.2.1], to ensure that the penalized residual sum of squares decreases at each iteration.
Convergence is declared when the orthogonality convergence criterion [@bateswatts88:_nonlin Sect. 2.2.3] is below some pre-specified tolerance.

The *Laplace approximation* to the deviance is
$$
  d(\mathbf{\theta},\mathbf{\beta},\sigma|\mathbf{y}_{\text{obs}})\approx
  n\log(2\pi\sigma^2)+2\log|\mathbf{L}_{\theta,\beta}|+
  \frac{r^2_{\theta,\beta}}{\sigma^2},
$$ {#eq:Laplace}
where the Cholesky factor, $\mathbf{L}_{\theta,\beta}$, and the penalized residual sum of squares, $r^2_{\theta,\beta}$, are both evaluated at the conditional mode, $\tilde{\mathbf{u}}$.
The Cholesky factor depends on $\mathbf{\theta}$, $\mathbf{\beta}$ and $\mathbf{u}$ for these models but typically the dependence on $\mathbf{\beta}$ and $\mathbf{u}$ is weak.

## Chapter summary {#sec:lmmsummary}

The definitions and the computational results for maximum likelihood estimation of the parameters in linear mixed models were summarized in @sec:defnLMM.
A key computation is evaluation of the blocked Cholesky factor, $\Lambda_\theta$, satisfying
@eq:blockedCholeskyP,
$$
\mathbf{L}_\theta\mathbf{L}_\theta'
   = \left(\Lambda_\theta'\mathbf{Z}'\mathbf{Z}\Lambda_\theta
   +\mathbf{I}_q\right) ,
$$

An extended decomposition @eq:bigdecomp provides the $q\times p$ matrix $\mathbf{R}_{ZX}$ and the $p\times p$ upper triangular $\mathbf{R}_X$ that are used to determine the conditional mode $\tilde{\mathbf{u}}_\theta$, the conditional estimate $\widehat{\mathbf{\beta}}_\theta$, and the minimum penalized residual sum of squares, $r^2_\theta$, from which the profiled deviance
$$
\tilde{d}(\mathbf{\theta}|\mathbf{y}_{\text{obs}})
  =2\log|\mathbf{L}_\theta|+n\left[1 +
    \log\left(\frac{2 \pi r^2_\theta}{n}\right)\right]
$$
or the profiled REML criterion
$$
\tilde{d}_R(\mathbf{\theta}|\mathbf{y}_{\text{obs}})=
  2\log\left(|\mathbf{L}_\theta||\mathbf{R}_X|\right)+(n-p)
  \left[1+\log\left(\frac{2\pi r^2_\theta}{n-p}\right)\right]
$$
can be evaluated and optimized (minimized) with respect to $\mathbf{\theta}$.
