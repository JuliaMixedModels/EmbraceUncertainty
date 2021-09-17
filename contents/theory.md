# Computational Methods for Mixed Models {#sec:computational}

In this chapter we describe some of the details of the computational methods for fitting linear mixed models, as implemented in the package, and the theoretical development behind these methods.
We also provide the basis for later generalizations to models for non-Gaussian responses and to models in which the relationship between the conditional mean, $\mathbf{\mu}$, and the linear predictor, $\mathbf{\gamma}=\mathbf{X}\mathbf{\beta}+\mathbf{Z}\mathbf{b}= \mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\mathbf{\beta}$, is a nonlinear relationship.

This material is directed at those readers who wish to follow the theory and methodology of linear mixed models and how both can be extended to other forms of mixed models.
Readers who are less interested in the "how" and the "why" of fitting mixed models than in the results themselves should not feel obligated to master these details.

We begin by reviewing the definition of linear mixed-effects models and some of the basics of the computational methods, as given in .

## Definitions and Basic Results {#sec:defnLMM}

As described in , a linear mixed-effects model is based on two vector-valued random variables: the $q$-dimensional vector of random effects, $\mathcal{B}$, and the $n$-dimensional response vector, $\mathcal{Y}$. Equation ([\[eq:LMMdist\]](#eq:LMMdist){reference-type="ref" reference="eq:LMMdist"}) defines the unconditional distribution of $\mathcal{B}$ and the conditional distribution of $\mathcal{Y}$, given $\mathcal{B}=\mathbf{b}$, as multivariate Gaussian distributions of the form
$$
\begin{aligned}
  (\mathcal{Y}|\mathcal{B}=\mathbf{b})&\sim\mathcal{N}(\mathbf{X}\mathbf{\beta}+\mathbf{Z}\mathbf{b},\sigma^2\mathbf{I})\\
  \mathcal{B}&\sim\mathcal{N}(\mathbf{0},\Sigma_\theta) .
\end{aligned}
$$

The $q\times q$, symmetric, variance-covariance matrix, $\mathrm{Var}(\mathcal{B})=\Sigma_\theta$, depends on the *variance-component parameter vector*, $\mathbf{\theta}$, and is *positive semidefinite*, which means that
$$\label{eq:posSemiDef} 
  \mathbf{b}'\Sigma_\theta\mathbf{b}\ge0,\quad\forall\,\mathbf{b}\ne\mathbf{0} .
$$
(The symbol $\forall$ denotes "for all".)
The fact that $\Sigma_\theta$ is positive semidefinite does not guarantee that $\Sigma_\theta^{-1}$ exists.
We would need a stronger property, $\mathbf{b}'\Sigma_\theta\mathbf{b}>0,\,\forall\,\mathbf{b}\ne\mathbf{0}$, called *positive definiteness*, to ensure that $\Sigma_\theta^{-1}$ exists.

Many computational formulas for linear mixed models are written in terms of the *precision matrix*, $\Sigma_\theta^{-1}$.
Such formulas will become unstable as $\Sigma_\theta$ approaches singularity.
And it can do so.
It is a fact that singular (i.e. non-invertible) $\Sigma_\theta$ can and do occur in practice, as we have seen in some of the examples in earlier chapters.
Moreover, during the course of the numerical optimization by which the parameter estimates are determined, it is frequently the case that the deviance or the REML criterion will need to be evaluated at values of $\mathbf{\theta}$ that produce a singular $\Sigma_\theta$.
Because of this we will take care to use computational methods that can be applied even when $\Sigma_\theta$ is singular and are stable as $\Sigma_\theta$ approaches singularity.

As defined in ([\[eq:relcovfac\]](#eq:relcovfac){reference-type="ref" reference="eq:relcovfac"}) a relative covariance factor, $\Lambda_\theta$, is any matrix that satisfies
$$
\Sigma_\theta=\sigma^2\Lambda_\theta\Lambda_\theta' .
$$
According to this definition, $\Sigma$ depends on both $\sigma$ and $\theta$, and we should write it as $\Sigma_{\sigma,\theta}$.
However, we will blur that distinction and continue to write $\text{Var}(\mathcal{B})=\Sigma_\theta$.
Another technicality is that the *common scale parameter*, $\sigma$, can, in theory, be zero.
We will show that in practice the only way for its estimate, $\widehat{\sigma}$, to be zero is for the fitted values from the fixed-effects only, $\mathbf{X}\widehat{\mathbf{\beta}}$, to be exactly equal to the observed data.
This occurs only with data that have been (incorrectly) simulated without error.
In practice we can safely assume that $\sigma>0$.
However, $\Lambda_\theta$, like $\Sigma_\theta$, can be singular.

Our computational methods are based on $\Lambda_\theta$ and do not require evaluation of $\Sigma_\theta$.
In fact, $\Sigma_\theta$ is explicitly evaluated only at the converged parameter estimates.

The spherical random effects, $\mathcal{U}\sim\mathcal{N}(\mathbf{0},\sigma^2\mathbf{I}_q)$, determine $\mathcal{B}$ as
$$\label{eq:sphericalRE}
  \mathcal{B}=\Lambda_\theta\mathcal{U} .
$$
Although it may seem more intuitive to write $\mathcal{U}$ as a linear transformation of $\mathcal{B}$, we cannot do that when $\Lambda_\theta$ is singular, which is why ([\[eq:sphericalRE\]](#eq:sphericalRE){reference-type="ref" reference="eq:sphericalRE"}) is in the form shown.

We can easily verify that ([\[eq:sphericalRE\]](#eq:sphericalRE){reference-type="ref" reference="eq:sphericalRE"}) provides the desired distribution for $\mathcal{B}$.
As a linear transformation of a multivariate Gaussian random variable, $\mathcal{B}$ will also be multivariate Gaussian.
Its mean and variance-covariance matrix are straightforward to evaluate,
$$\label{eq:EB}
  \mathrm{E}[\mathcal{B}] = \Lambda_\theta\mathrm{E}[\mathcal{U}]=\Lambda_\theta\mathbf{0}=\mathbf{0}
$$
and
$$\mathrm{Var}(\mathcal{B})
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
$$\label{eq:condYgivenU}
  \begin{aligned}
  (\mathcal{Y}|\mathcal{U}=\mathbf{u})&\sim\mathcal{N}(\mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\beta,\sigma^2\mathbf{I}_n)\\
  \mathcal{U}&\sim\mathcal{N}(\mathbf{0},\sigma^2\mathbf{I}_q) .
  \end{aligned}
$$

To allow for extensions to other types of mixed models we distinguish between the *linear predictor* 
$$\label{eq:linearpred}
  \mathbf{\gamma} = \mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\beta
$$
and the *conditional mean* of $\mathcal{Y}$, given $\mathcal{U}=\mathbf{u}$, which is
$$\label{eq:conditionalMean}
  \mathbf{\mu} = \mathrm{E}\left[\mathcal{Y}|\mathcal{U}=\mathbf{u}\right] .
$$
For a linear mixed model $\mathbf{\mu}=\mathbf{\gamma}$.
In other forms of mixed models the conditional mean, $\mathbf{\mu}$, can be a nonlinear function of the linear predictor, $\mathbf{\gamma}$.
For some models the dimension of $\mathbf{\gamma}$ is a multiple of $n$, the dimension of $\mathbf{\mu}$ and $\mathbf{y}$, but for a linear mixed model the dimension of $\mathbf{\gamma}$ must be $n$.
Hence, the model matrix $\mathbf{Z}$ must be $n\times q$ and $\mathbf{X}$ must be $n\times p$.

## The Conditional Distribution $(\mathcal{U}|\mathcal{Y}=\mathbf{y})$ {#sec:conddistUgivenY}

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
$$\label{eq:conddenUgivenY}
  f_{\mathcal{U}|\mathcal{Y}}(\mathbf{u}|\mathbf{y}_{\text{obs}})=\frac
  {h(\mathbf{u})}{\int_{\mathbb{R}^q}h(\mathbf{u})\,d\mathbf{u}}  .
$$
In some theoretical developments the normalizing constant, which is the integral
in the denominator of an expression like ([\[eq:conddenUgivenY\]](#eq:conddenUgivenY){reference-type="ref" reference="eq:conddenUgivenY"}), is not of interest.
Here it is of interest because the normalizing constant is exactly the likelihood that we wish to evaluate,
$$\label{eq:LMMlikelihood}
  L(\mathbf{\theta},\mathbf{\beta},\sigma|\mathbf{y}_{\text{obs}}) =
  \int_{\mathbb{R}^q}h(\mathbf{u})\,d\mathbf{u} .
$$

For a linear mixed model, where all the distributions of interest are
multivariate Gaussian and the conditional mean, $\mathbf{\mu}$, is a linear
function of both $\mathbf{u}$ and $\mathbf{\beta}$, the distinction between
evaluating the joint density at $\mathbf{y}_{\text{obs}}$ to produce
$h(\mathbf{u})$ then integrating with respect to $\mathbf{u}$, as opposed to
first integrating the joint density then evaluating at $\mathbf{y}_{\text{obs}}$, is not terribly important.
For other mixed models this distinction can be important.
In particular, generalized linear mixed models, described in Chap. [\[chap:GLMMbinomial\]](#chap:GLMMbinomial){reference-type="ref" reference="chap:GLMMbinomial"}, are often used to model a discrete response, such as a binary response or a count, leading to a joint distribution for $\mathcal{Y}$ and $\mathcal{U}$ that is discrete with respect to one variable, $\mathbf{y}$, and continuous with respect to the other, $\mathbf{u}$.
In such cases there isn't a joint density for $\mathcal{Y}$ and $\mathcal{U}$. The necessary distribution theory for general $\mathbf{y}$ and $\mathbf{u}$ is well-defined but somewhat awkward to describe.
It is much easier to realize that we are only interested in the observed response vector, $\mathbf{y}_{\text{obs}}$, not some arbitrary value of $\mathbf{y}$, so we can concentrate on the conditional distribution of $\mathcal{U}$ given $\mathcal{Y}=\mathbf{y}_{\text{obs}}$. For all the mixed models we will consider, the conditional distribution, $(\mathcal{U}|\mathcal{Y}=\mathbf{y}_{\text{obs}})$, is continuous and both the conditional density, $f_{\mathcal{U}|\mathcal{Y}}(\mathbf{u}|\mathbf{y}_{\text{obs}})$, and its unnormalized form, $h(\mathbf{u})$, are well-defined.

## Integrating $h(\mathbf{u})$ in the Linear Mixed Model {#sec:IntegratingH}

The integral defining the likelihood in ([\[eq:LMMlikelihood\]](#eq:LMMlikelihood){reference-type="ref" reference="eq:LMMlikelihood"}) has a closed form in the case of a linear mixed model but not for some of the more general forms of mixed models.
To motivate methods for approximating the likelihood in more general situations, we describe in some detail how the integral can be evaluated using the sparse Cholesky factor, $\mathbf{L}_\theta$, and the conditional mode, $$\label{eq:condMode}
  \tilde{\mathbf{u}}=\arg\max_{\mathbf{u}} f_{\mathcal{U}|\mathcal{Y}}(\mathbf{u}|\mathbf{y}_{\text{obs}})=
  \arg\max_{\mathbf{u}} h(\mathbf{u}) = \arg\max_{\mathbf{u}}
  f_{\mathcal{Y}|\mathcal{U}}(\mathbf{y}_{\text{obs}}|\mathbf{u})\,f_{\mathcal{U}}(\mathbf{u}).
$$
The notation $\arg\max_{\mathbf{u}}$ means that $\tilde{\mathbf{u}}$ is the value of $\mathbf{u}$ that maximizes the expression that follows.

In general, the *mode* of a continuous distribution is the value of the random variable that maximizes the density.
The value $\tilde{\mathbf{u}}$ is called the conditional mode of $\mathbf{u}$, given $\mathcal{Y}=\mathbf{y}_{\text{obs}}$, because $\tilde{\mathbf{u}}$ maximizes the conditional density of $\mathcal{U}$ given $\mathcal{Y}=\mathbf{y}_{\text{obs}}$.
The location of the maximum can be determined by maximizing the unnormalized
conditional density because $h(\mathbf{u})$ is just a constant multiple of
$f_{\mathcal{U}|\mathcal{Y}}(\mathbf{u}|\mathbf{y}_{\text{obs}})$.
The last part of ([\[eq:condMode\]](#eq:condMode){reference-type="ref" reference="eq:condMode"}) is simply a re-expression of $h(\mathbf{u})$ as the product of $f_{\mathcal{Y}|\mathcal{U}}(\mathbf{y}_{\text{obs}}|\mathbf{u})$ and $f_{\mathcal{U}}(\mathbf{u})$.
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
$$\label{eq:hudef}
  h(\mathbf{u})=\frac{1}{\left(2\pi\sigma^2\right)^{(n+q)/2}}
  \exp\left(-\frac{\left\|\mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+\|\mathbf{u}\|^2}{2\sigma^2}\right) .
$$
On the deviance scale we have
$$\label{eq:devh}
  -2\log\left(h(\mathbf{u})\right)=(n+q)\log(2\pi\sigma^2)
  +\frac{\left\|\mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+\|\mathbf{u}\|^2}{\sigma^2} .
$$
Because ([\[eq:devh\]](#eq:devh){reference-type="ref" reference="eq:devh"})
describes the negative log density, $\tilde{\mathbf{u}}$ will be the value
of $\mathbf{u}$ that minimizes the expression on the right hand side of
([\[eq:devh\]](#eq:devh){reference-type="ref" reference="eq:devh"}).

The only part of the right hand side of ([\[eq:devh\]](#eq:devh){reference-type="ref" reference="eq:devh"}) that depends on $\mathbf{u}$ is the numerator of the second term.
Thus 
$$\label{eq:PLSsol}
  \tilde{\mathbf{u}}=\arg\min_{\mathbf{u}} \left\|
    \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+
  \|\mathbf{u}\|^2.
$$
The expression to be minimized, called the *objective function*, is described as a *penalized residual sum of squares* (PRSS) and the minimizer, $\tilde{\mathbf{u}}$, is called the *penalized least squares* (PLS) solution.
They are given these names because the first term in the objective, $\left\| \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2$, is a sum of squared residuals, and the second term, $\|\mathbf{u}\|^2$, is a penalty on the length, $\|\mathbf{u}\|$, of $\mathbf{u}$.
Larger values of $\mathbf{u}$ (in the sense of greater lengths as vectors) incur a higher penalty.

The PRSS criterion determining the conditional mode balances fidelity to the observed data (i.e. producing a small residual sum of squares) against simplicity of the model (small $\|\mathbf{u}\|$).
We refer to this type of criterion as a smoothing objective, in the sense that it seeks to smooth out the fitted response by reducing model complexity while
still retaining reasonable fidelity to the observed data.

For the purpose of evaluating the likelihood we will regard the PRSS criterion as a function of the parameters, given the data, and write its minimum value as $$\label{eq:r2thetabeta}
  r^2_{\theta,\beta}=\min_{\mathbf{u}} \left\|
    \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+ \|\mathbf{u}\|^2.
$$
Notice that $\mathbf{\beta}$ only enters the right hand side of ([\[eq:r2thetabeta\]](#eq:r2thetabeta){reference-type="ref" reference="eq:r2thetabeta"}) through the linear predictor expression.
We will see that $\tilde{\mathbf{u}}$ can be determined by a direct (i.e. non-iterative) calculation and, in fact, we can minimize the PRSS
criterion with respect to $\mathbf{u}$ and $\mathbf{\beta}$ simultaneously without iterating.
We write this minimum value as
$$\label{eq:r2theta}
  r^2_\theta=\min_{\mathbf{u},\mathbf{\beta}} \left\|
    \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+ \|\mathbf{u}\|^2.
$$
The value of $\mathbf{\beta}$ at the minimum is called the conditional estimate of $\mathbf{\beta}$ given $\mathbf{\theta}$, written $\widehat{\mathbf{\beta}}_\theta$.

## Determining the PLS Solutions, $\tilde{\mathbf{u}}$ and $\widehat{\mathbf{\beta}}_\theta$ {#sec:PLSsol}

One way of expressing a penalized least squares problem like ([\[eq:r2thetabeta\]](#eq:r2thetabeta){reference-type="ref" reference="eq:r2thetabeta"}) is by incorporating the penalty as "pseudo-data" in an ordinary least squares problem.
We extend the "response vector", which is $\mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}$ when we minimize with respect to $\mathbf{u}$ only, with $q$ responses that are 0 and we extend the predictor expression, $\mathbf{Z}\Lambda_\theta\mathbf{u}$ with $\mathbf{I}_q\mathbf{u}$.
Writing this as a least squares problem produces 
$$\label{eq:PLSLMM}
  \tilde{\mathbf{u}}=\arg\min_{\mathbf{u}}\left\|
    \begin{bmatrix}
      \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}\\
      \mathbf{0}
    \end{bmatrix} -
    \begin{bmatrix}
      \mathbf{Z}\Lambda_\theta\\
      \mathbf{I}_q
    \end{bmatrix}\mathbf{u}\right\|^2
$$
with a solution that satisfies
$$\label{eq:LMMPLSsol}
  \left(\Lambda_\theta'\mathbf{Z}'\mathbf{Z}\Lambda_\theta+\mathbf{I}_q\right)\tilde{\mathbf{u}}
  =\Lambda_\theta'\mathbf{Z}'\left(\mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}\right).
$$

To evaluate $\tilde{\mathbf{u}}$ we form the *sparse Cholesky factor*, $\mathbf{L}_\theta$, which is a lower triangular $q\times q$ matrix that satisfies $$\label{eq:sparseCholesky}
  \mathbf{L}_\theta\mathbf{L}_\theta'=
  \Lambda_\theta'\mathbf{Z}'\mathbf{Z}\Lambda_\theta+\mathbf{I}_q .
$$
The actual evaluation of the sparse Cholesky factor, $\mathbf{L}_\theta$, often incorporates a *fill-reducing permutation*, which we describe next.

### The Fill-reducing Permutation, $\mathbf{P}$ {#sec:fill-reducingP}

In earlier chapters we have seen that often the random effects vector is
re-ordered before $\mathbf{L}_\theta$ is created. The re-ordering or
permutation of the elements of $\mathbf{u}$ and, correspondingly, the
columns of the model matrix, $\mathbf{Z}\Lambda_\theta$, does not affect the
theory of linear mixed models but can have a profound effect on the time
and storage required to evaluate $\mathbf{L}_\theta$ in large problems. We
write the effect of the permutation as multiplication by a $q\times q$
*permutation matrix*, $\mathbf{P}$, although in practice we apply the
permutation without ever constructing $\mathbf{P}$.
That is, the matrix $\mathbf{P}$ is a notational convenience only.

The matrix $\mathbf{P}$ consists of permuted columns of the identity matrix,
$\mathbf{I}_q$, and it is easy to establish that the inverse permutation
corresponds to multiplication by $\mathbf{P}'$. Because multiplication
by $\mathbf{P}$ or by $\mathbf{P}'$ simply re-orders the components of a
vector, the length of the vector is unchanged. Thus,
$$\label{eq:orthogonalP}
  \|\mathbf{P}\mathbf{u}\|^2= \|\mathbf{u}\|^2 = \|\mathbf{P}'\mathbf{u}\|^2$$ and we
can express the penalty in
([\[eq:r2theta\]](#eq:r2theta){reference-type="ref"
reference="eq:r2theta"}) in any of these three forms. The properties of
$\mathbf{P}$ that it preserves lengths of vectors and that its transpose is
its inverse are summarized by stating that $\mathbf{P}$ is an *orthogonal
matrix*.

The permutation represented by $\mathbf{P}$ is determined from the structure
of $\Lambda_\theta'\mathbf{Z}'\mathbf{Z}\Lambda_\theta+\mathbf{I}_q$ for some initial value of $\mathbf{\theta}$. The
particular value of $\mathbf{\theta}$ does not affect the result because the
permutation depends only the positions of the non-zeros, not the
numerical values at these positions.

Taking into account the permutation, the sparse Cholesky factor, $\mathbf{L}_\theta$, is defined to be the sparse, lower triangular, $q\times q$ matrix with positive diagonal elements satisfying
$$
\label{eq:sparseCholeskyP}
  \mathbf{L}_\theta\mathbf{L}_\theta'
   = \mathbf{P}\left(\Lambda_\theta'\mathbf{Z}'\mathbf{Z}\Lambda_\theta
   +\mathbf{I}_q\right)\mathbf{P}'.
$$
(Problems [\[pr:th:posdef\]](#pr:th:posdef){reference-type="ref"
reference="pr:th:posdef"} and
[\[pr:th:posdiag\]](#pr:th:posdiag){reference-type="ref"
reference="pr:th:posdiag"} outline the steps in showing that we can
require the diagonal elements of $\Lambda_\theta$ to be positive, not
just non-negative.)
Problems [\[pr:th:posdef\]](#pr:th:posdef){reference-type="ref"
reference="pr:th:posdef"} and
[\[pr:th:posdiag\]](#pr:th:posdiag){reference-type="ref"
reference="pr:th:posdiag"} indicate why we can require this. Because the
diagonal elements of $\Lambda_\theta$ are positive, its determinant,
$|\Lambda_\theta|$, which, for a triangular matrix such as
$\Lambda_\theta$, is simply the product of its diagonal elements, is
also positive.

Many sparse matrix methods, including the sparse Cholesky decomposition, are performed in two stages: the *symbolic phase* in which the locations of the non-zeros in the result are determined and the *numeric phase* in which the numeric values at these positions are evaluated.
The symbolic phase for the decomposition
([\[eq:sparseCholeskyP\]](#eq:sparseCholeskyP){reference-type="ref"
reference="eq:sparseCholeskyP"}), which includes determining the
permutation, $\mathbf{P}$, need only be done once. Evaluation of $\mathbf{L}_\theta$ for subsequent values of $\mathbf{\theta}$ requires only the
numeric phase, which typically is much faster than the symbolic phase.

The permutation, $\mathbf{P}$, serves two purposes. The first, and more
important purpose, is to reduce the number of non-zeros in the factor,
$\mathbf{L}_\theta$. The factor is potentially non-zero at every non-zero
location in the lower triangle of the matrix being decomposed. However,
as we saw in
Fig. [\[fig:fm03LambdaLimage\]](#fig:fm03LambdaLimage){reference-type="ref"
reference="fig:fm03LambdaLimage"} of , there may be positions in the
factor that get filled-in even though they are known to be zero in the
matrix being decomposed. The *fill-reducing permutation* is chosen
according to certain heuristics to reduce the amount of fill-in. We use
the approximate minimal degree (AMD) method described in @Davis:1996.
After the fill-reducing permutation is determined, a "post-ordering" is
applied. This has the effect of concentrating the non-zeros near the
diagonal of the factor. See @davis06:csparse_book for details.

The pseudo-data representation of the PLS problem,
([\[eq:PLSLMM\]](#eq:PLSLMM){reference-type="ref"
reference="eq:PLSLMM"}), becomes $$\label{eq:PLSLMMP}
  \tilde{\mathbf{u}}=\arg\min_{\mathbf{u}}\left\|
    \begin{bmatrix}
      \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}\\
      \mathbf{0}
    \end{bmatrix} -
    \begin{bmatrix}
      \mathbf{Z}\Lambda_\theta\mathbf{P}'\\
      \mathbf{P}'
    \end{bmatrix}\mathbf{P}\mathbf{u}\right\|^2$$ and the system of linear
equations satisfied by $\tilde{\mathbf{u}}$ is $$\label{eq:LMMPLSsolP}
  \mathbf{L}_\theta\mathbf{L}_\theta'\mathbf{P}\tilde{\mathbf{u}}=
  \mathbf{P}\left(\Lambda_\theta'\mathbf{Z}'\mathbf{Z}\Lambda_\theta+\mathbf{I}_q\right)\mathbf{P}'\mathbf{P}\tilde{\mathbf{u}}
  =\mathbf{P}\Lambda_\theta'\mathbf{Z}'\left(\mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}\right) .$$

Obtaining the Cholesky factor, $\mathbf{L}_\theta$, may not seem to be great
progress toward determining $\tilde{\mathbf{u}}$ because we still must solve
([\[eq:LMMPLSsolP\]](#eq:LMMPLSsolP){reference-type="ref"
reference="eq:LMMPLSsolP"}) for $\tilde{\mathbf{u}}$. However, it is the key
to the computational methods in the package. The ability to evaluate
$\mathbf{L}_\theta$ rapidly for many different values of $\mathbf{\theta}$ is
what makes the computational methods in feasible, even when applied to
very large data sets with complex structure. Once we evaluate
$\mathbf{L}_\theta$ it is straightforward to solve
([\[eq:LMMPLSsolP\]](#eq:LMMPLSsolP){reference-type="ref"
reference="eq:LMMPLSsolP"}) for $\tilde{\mathbf{u}}$ because $\mathbf{L}_\theta$
is triangular.

In we will describe the steps in determining this solution but first we
will show that the solution, $\tilde{\mathbf{u}}$, and the value of the
objective at the solution, $r^2_{\theta,\beta}$, do allow us to evaluate
the deviance.

### The Value of the Deviance and Profiled Deviance {#sec:solvingtildeu}

After evaluating $\mathbf{L}_\theta$ and using that to solve for
$\tilde{\mathbf{u}}$, which also produces $r^2_{\beta,\theta}$, we can write
the PRSS for a general $\mathbf{u}$ as $$\label{eq:PRSSwithL}
\left\|
    \mathbf{y}_{\text{obs}}-\mathbf{X}\mathbf{\beta}-\mathbf{Z}\Lambda_\theta\mathbf{u}\right\|^2+ \|\mathbf{u}\|^2=
  r^2_{\theta,\beta}+\|\mathbf{L}_\theta'(\mathbf{u}-\tilde{\mathbf{u}})\|^2$$
which finally allows us to evaluate the likelihood. We plug the right
hand side of ([\[eq:PRSSwithL\]](#eq:PRSSwithL){reference-type="ref"
reference="eq:PRSSwithL"}) into the definition of $h(\mathbf{u})$ and apply
the change of variable $$\label{eq:changeVar}
  \mathbf{Z}=\frac{\mathbf{L}_\theta'(\mathbf{u}-\tilde{\mathbf{u}})}{\sigma} .$$
The determinant of the Jacobian of this transformation,
$$\label{eq:tranJac}
  \left|\frac{d\mathbf{Z}}{d\mathbf{u}}\right|=
  \left|\frac{\mathbf{L}_\theta'}{\sigma}\right|=
  \frac{|\mathbf{L}_\theta|}{\sigma^q}$$ is required for the change of
variable in the integral. We use the letter $\mathbf{Z}$ for the transformed
value because we will rearrange the integral to have the form of the
integral of the density of the standard multivariate normal
distribution. That is, we will use the result
$$\label{eq:stdnormint}
   \int_{\mathbb{R}^q}\frac{e^{-\|\mathbf{Z}\|^2/2}}{(2\pi)^{q/2}}\,d\mathbf{z} = 1.
$$

Putting all these pieces together gives $$\label{eq:hintegral}
  \begin{aligned}
    L(\theta,\beta,\sigma)&=\int_{\mathbb{R}^q}h(\mathbf{u})\,d\mathbf{u}\\
    &=\int_{\mathbb{R}^q}\frac{1}{(2\pi\sigma^2)^{(n+q)/2}}
    \exp\left(-\frac{r^2_{\theta,\beta}+\|\mathbf{L}_\theta'(\mathbf{u}-\tilde{\mathbf{u}})\|^2}{2\sigma^2}\right)d\mathbf{u}\\
    &=\frac{\exp\left(-\frac{r^2_{\theta,\beta}}{2\sigma^2}\right)}
    {(2\pi\sigma^2)^{n/2}}\int_{\mathbb{R}^q}\frac{1}{(2\pi)^{q/2}}
    \exp\left(-\frac{\|\mathbf{L}_\theta'(\mathbf{u}-\tilde{\mathbf{u}})\|^2}{2\sigma^2}\right)
    \frac{|\mathbf{L}_\theta|}{|\mathbf{L}_\theta|}
    \frac{d\mathbf{u}}{\sigma^q}\\
    &=\frac{\exp\left(-\frac{r^2_{\theta,\beta}}{2\sigma^2}\right)}
    {(2\pi\sigma^2)^{n/2}|\mathbf{L}_\theta|}
    \int_{\mathbb{R}^q}\frac{e^{-\|\mathbf{z}\|^2/2}}{(2\pi)^{q/2}}\,d\mathbf{Z}\\
    &=\frac{\exp\left(-\frac{r^2_{\theta,\beta}}{2\sigma^2}\right)}
    {(2\pi\sigma^2)^{n/2}|\mathbf{L}_\theta|} .
  \end{aligned}$$

The deviance can now be expressed as
$$d(\mathbf{\theta},\mathbf{\beta},\sigma|\mathbf{y}_{\text{obs}})=
  -2\log\left(L(\mathbf{\theta},\mathbf{\beta},\sigma|\mathbf{y}_{\text{obs}})\right)
  =n\log(2\pi\sigma^2)+2\log|\mathbf{L}_\theta|+
  \frac{r^2_{\beta,\theta}}{\sigma^2},$$ as stated in
([\[eq:LMMdeviance\]](#eq:LMMdeviance){reference-type="ref"
reference="eq:LMMdeviance"}). The maximum likelihood estimates of the
parameters are those that minimize this deviance.

Equation ([\[eq:LMMdeviance\]](#eq:LMMdeviance){reference-type="ref"
reference="eq:LMMdeviance"}) is a remarkably compact expression,
considering that the class of models to which it applies is very large
indeed. However, we can do better than this if we notice that
$\mathbf{\beta}$ affects
([\[eq:LMMdeviance\]](#eq:LMMdeviance){reference-type="ref"
reference="eq:LMMdeviance"}) only through $r^2_{\beta,\theta}$, and, for
any value of $\mathbf{\theta}$, minimizing this expression with respect to
$\mathbf{\beta}$ is just an extension of the penalized least squares problem.
Let $\widehat{\mathbf{\beta}}_\theta$ be the value of $\mathbf{\beta}$ that
minimizes the PRSS simultaneously with respect to $\mathbf{\beta}$ and
$\mathbf{u}$ and let $r^2_\theta$ be the PRSS at these minimizing values.
If, in addition, we set $\widehat{\sigma^2}_\theta=r^2_\theta/n$, which
is the value of $\sigma^2$ that minimizes the deviance for a given value
of $r^2_\theta$, then the *profiled deviance*, which is a function of
$\mathbf{\theta}$ only, becomes $$\label{eq:LMMprofdeviance}
  \tilde{d}(\mathbf{\theta}|\mathbf{y}_{\text{obs}})
  =2\log|\mathbf{L}_\theta|+n\left[1 +
    \log\left(\frac{2 \pi r^2_\theta}{n}\right)\right].$$

Numerical optimization (minimization) of $\tilde{d}(\mathbf{\theta}|\mathbf{y}_{\text{obs}})$ with respect to $\mathbf{\theta}$ determines the MLE,
$\widehat{\mathbf{\theta}}$. The MLEs for the other parameters,
$\widehat{\mathbf{\beta}}$ and $\widehat{\sigma}$, are the corresponding
conditional estimates evaluated at $\widehat{\mathbf{\theta}}$.

### Determining $r^2_\theta$ and $\hat{\mathbf{\beta}}_\theta$ {#sec:betahat}

To determine $\tilde{\mathbf{u}}$ and $\widehat{\mathbf{\beta}}_\theta$
simultaneously we rearrange the terms in
([\[eq:PLSLMMP\]](#eq:PLSLMMP){reference-type="ref"
reference="eq:PLSLMMP"}) as $$\label{eq:PLSLMM1}
  \begin{bmatrix}
    \tilde{\mathbf{u}}\\
    \widehat{\mathbf{\beta}}_\theta
  \end{bmatrix}
  =\arg\min_{\mathbf{u},\mathbf{\beta}}
  \left\|
    \begin{bmatrix}
      \mathbf{y}_{\text{obs}}\\
      \mathbf{0}
    \end{bmatrix} -
    \begin{bmatrix}
      \mathbf{Z}\Lambda_\theta\mathbf{P}' & \mathbf{X}\\
      \mathbf{P}' &\mathbf{0}
    \end{bmatrix}
    \begin{bmatrix}
      \mathbf{P}\mathbf{u}\\
      \mathbf{\beta}
    \end{bmatrix}
  \right\|^2 .$$ The PLS values, $\tilde{\mathbf{u}}$ and
$\widehat{\mathbf{\beta}}_\theta$, are the solutions to $$\label{eq:bigPLS}
  \begin{bmatrix}
    \mathbf{P}\left(\Lambda_\theta'\mathbf{Z}'\mathbf{Z}\Lambda_\theta+\mathbf{I}_q\right)\mathbf{P}' &
    \mathbf{P}\Lambda_\theta'\mathbf{Z}'\mathbf{X}\\
    \mathbf{X}'\mathbf{Z}\Lambda_\theta\mathbf{P}' & \mathbf{X}'\mathbf{X}
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{P}\tilde{\mathbf{u}}\\
    \widehat{\mathbf{\beta}}_\theta
  \end{bmatrix}=
  \begin{bmatrix}
    \mathbf{P}\Lambda_\theta'\mathbf{Z}'\mathbf{y}_{\text{obs}}\\
    \mathbf{X}'\mathbf{y}_{\text{obs}}
  \end{bmatrix}.$$ To evaluate these solutions we decompose the system
matrix as $$\label{eq:bigdecomp}
  \begin{bmatrix}
    \mathbf{P}\left(\Lambda_\theta'\mathbf{Z}'\vec
      Z\Lambda_\theta+\mathbf{I}_q\right)\mathbf{P}' &
    \mathbf{P}\Lambda_\theta'\mathbf{Z}'\mathbf{X}\\
    \mathbf{X}'\mathbf{Z}\Lambda_\theta\mathbf{P}' & \mathbf{X}'\mathbf{X}
  \end{bmatrix}
  =
  \begin{bmatrix}
    \mathbf{L}_\theta & \mathbf{0}\\
    \mathbf{R}_{ZX}' & \mathbf{R}_X'
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{L}_\theta' & \mathbf{R}_{ZX}\\
    \mathbf{0} & \mathbf{R}_X
  \end{bmatrix}$$ where, as before, $\mathbf{L}_\theta$, the sparse Cholesky
factor, is the sparse lower triangular $q\times q$ matrix satisfying
([\[eq:sparseCholeskyP\]](#eq:sparseCholeskyP){reference-type="ref"
reference="eq:sparseCholeskyP"}). The other two matrices in
([\[eq:bigdecomp\]](#eq:bigdecomp){reference-type="ref"
reference="eq:bigdecomp"}): $\mathbf{R}_{ZX}$, which is a general
$q\times p$ matrix, and $\mathbf{R}_X$, which is an upper triangular
$p\times p$ matrix, satisfy $$\label{eq:RZXdef}
  \mathbf{L}_\theta\mathbf{R}_{ZX}=\mathbf{P}\Lambda_\theta'\mathbf{Z}'\mathbf{X}$$
and $$\label{eq:RXdef}
  \mathbf{R}_X'\mathbf{R}_X=\mathbf{X}'\mathbf{X}-\mathbf{R}_{ZX}'\mathbf{R}_{ZX}.$$

Those familiar with standard ways of writing a Cholesky decomposition as
either $\mathbf{L}\mathbf{L}'$ or $\mathbf{R}'\mathbf{R}$ ($\mathbf{L}$ is the
factor as it appears on the left and $\mathbf{R}$ is as it appears on the
right) will notice a notational inconsistency in
([\[eq:bigdecomp\]](#eq:bigdecomp){reference-type="ref"
reference="eq:bigdecomp"}). One Cholesky factor is defined as the lower
triangular fractor on the left and the other is defined as the upper
triangular factor on the right. It happens that in the Cholesky factor
of a dense positive-definite matrix is returned as the right factor,
whereas the sparse Cholesky factor is returned as the left factor.

One other technical point that should be addressed is whether $\vec
X'\mathbf{X}-\mathbf{R}_{ZX}'\mathbf{R}_{ZX}$ is positive definite. In
theory, if $\mathbf{X}$ has full column rank, so that $\mathbf{X}'\mathbf{X}$
is positive definite, then the downdated matrix, $\mathbf{X}'\vec
X-\mathbf{R}_{ZX}'\mathbf{R}_{ZX}$, must also be positive definite (see
Prob. [\[pr:th:downdate\]](#pr:th:downdate){reference-type="ref"
reference="pr:th:downdate"}). In practice, the downdated matrix can
become computationally singular in ill-conditioned problems, in which
case an error is reported.

The extended decomposition
([\[eq:bigdecomp\]](#eq:bigdecomp){reference-type="ref"
reference="eq:bigdecomp"}) not only provides for the evaluation of the
profiled deviance function, $\tilde{d}(\mathbf{\theta})$,
([\[eq:LMMprofdeviance\]](#eq:LMMprofdeviance){reference-type="ref"
reference="eq:LMMprofdeviance"}) but also allows us to define and
evaluate the profiled REML criterion.

The REML Criterion {#sec:REML}
------------------

The so-called REML estimates of variance components are often preferred
to the maximum likelihood estimates. ("REML" can be considered to be an
acronym for "restricted" or "residual" maximum likelihood, although
neither term is completely accurate because these estimates do not
maximize a likelihood.) We can motivate the use of the REML criterion by
considering a linear regression model, $$\label{eq:20}
  \mathcal{Y}\sim\mathcal{N}(\mathbf{X}\mathbf{\beta},\sigma^2\mathbf{I}_n),$$ in which we
typically estimate $\sigma^2$ as $$\label{eq:21}
  \widehat{\sigma^2_R}=\frac{\|\mathbf{y}_{\text{obs}}-\mathbf{X}\widehat{\mathbf{\beta}}\|^2}{n-p}$$
even though the maximum likelihood estimate of $\sigma^2$ is
$$\label{eq:22}
  \widehat{\sigma^2_{L}}=\frac{\|\mathbf{y}_{\text{obs}}-\vec
    X\widehat{\mathbf{\beta}}\|^2}{n} .$$

The argument for preferring $\widehat{\sigma^2_R}$ to
$\widehat{\sigma^2_{L}}$ as an estimate of $\sigma^2$ is that the
numerator in both estimates is the sum of squared residuals at
$\widehat{\mathbf{\beta}}$ and, although the residual vector, $\vec
y_{\text{obs}}-\mathbf{X}\widehat{\mathbf{\beta}}$, is an $n$-dimensional vector,
it satisfies $p$ linearly independent constraints, $\vec
X'(\mathbf{y}_{\text{obs}}-\mathbf{X}\widehat{\mathbf{\beta}})=\mathbf{0}$. That is,
the residual at $\widehat{\mathbf{\beta}}$ is the projection of the observed
response vector, $\mathbf{y}_{\text{obs}}$, into an $(n-p)$-dimensional
linear subspace of the $n$-dimensional response space. The estimate
$\widehat{\sigma^2_R}$ takes into account the fact that $\sigma^2$ is
estimated from residuals that have only $n-p$ *degrees of freedom*.

Another argument often put forward for REML estimation is that
$\widehat{\sigma^2_R}$ is an *unbiased* estimate of $\sigma^2$, in the
sense that the expected value of the estimator is equal to the value of
the parameter. However, determining the expected value of an estimator
involves integrating with respect to the density of the estimator and we
have seen that densities of estimators of variances will be skewed,
often highly skewed. It is not clear why we should be interested in the
expected value of a highly skewed estimator. If we were to transform to
a more symmetric scale, such as the estimator of the standard deviation
or the estimator of the logarithm of the standard deviation, the REML
estimator would no longer be unbiased. Furthermore, this property of
unbiasedness of variance estimators does not generalize from the linear
regression model to linear mixed models. This is all to say that the
distinction between REML and ML estimates of variances and variance
components is probably less important than many people believe.

Nevertheless it is worthwhile seeing how the computational techniques
described in this chapter apply to the REML criterion because the REML
parameter estimates $\widehat{\mathbf{\theta}}_R$ and $\widehat{\sigma_R^2}$
for a linear mixed model have the property that they would specialize to
$\widehat{\sigma^2_R}$ from ([\[eq:21\]](#eq:21){reference-type="ref"
reference="eq:21"}) for a linear regression model, as seen in .

Although not usually derived in this way, the REML criterion (on the
deviance scale) can be expressed as $$\label{eq:23}
  d_R(\mathbf{\theta},\sigma|\mathbf{y}_{\text{obs}})=-2\log
  \int_{\mathbb{R}^p}L(\mathbf{\theta},\mathbf{\beta},\sigma|\mathbf{y}_{\text{obs}})\,d\mathbf{\beta} .$$
The REML estimates $\widehat{\mathbf{\theta}}_R$ and $\widehat{\sigma_R^2}$
minimize $d_R(\mathbf{\theta},\sigma|\mathbf{y}_{\text{obs}})$.

To evaluate this integral we form an expansion, similar to
([\[eq:PRSSwithL\]](#eq:PRSSwithL){reference-type="ref"
reference="eq:PRSSwithL"}), of $r^2_{\theta,\beta}$ about
$\widehat{\mathbf{\beta}}_\theta$ $$\label{eq:rsqbetathetaexp}
  r^2_{\theta,\beta}=r^2_\theta+\|\mathbf{R}_X(\mathbf{\beta}-\widehat{\mathbf{\beta}}_\theta)\|^2 .$$
In the same way that
([\[eq:PRSSwithL\]](#eq:PRSSwithL){reference-type="ref"
reference="eq:PRSSwithL"}) was used to simplify the integral in
([\[eq:hintegral\]](#eq:hintegral){reference-type="ref"
reference="eq:hintegral"}), we can derive $$\label{eq:betaintegral}
  \int_{\mathbb{R}^p}\frac{\exp\left(-\frac{r^2_{\theta,\beta}}{2\sigma^2}\right)}
  {(2\pi\sigma^2)^{n/2}|\mathbf{L}_\theta|} \,d\mathbf{\beta}=
  \frac{\exp\left(-\frac{r^2_\theta}{2\sigma^2}\right)}
  {(2\pi\sigma^2)^{(n-p)/2}|\mathbf{L}_\theta||\mathbf{R}_X|}$$ corresponding to
a REML criterion on the deviance scale of $$\label{eq:REMLdev}
  d_R(\mathbf{\theta},\sigma|\mathbf{y}_{\text{obs}})=(n-p)\log(2\pi\sigma^2)+
  2\log\left(|\mathbf{L}_\theta||\mathbf{R}_X|\right)+\frac{r^2_\theta}{\sigma^2} .$$
Plugging in the conditional REML estimate,
$\widehat{\sigma^2}_R=r^2_\theta/(n-p)$, provides the profiled REML
criterion $$\label{eq:24}
  \tilde{d}_R(\mathbf{\theta}|\mathbf{y}_{\text{obs}})=
  2\log\left(|\mathbf{L}_\theta||\mathbf{R}_X|\right)+(n-p)
  \left[1+\log\left(\frac{2\pi r^2_\theta}{n-p}\right)\right].$$

The REML estimate of $\mathbf{\theta}$ is $$\label{eq:31}
  \widehat{\mathbf{\theta}}_R=\arg\min_{\mathbf{\theta}}\tilde{d}_R(\mathbf{\theta}|\mathbf{y}_{\text{obs}}) ,$$
and the REML estimate of $\sigma^2$ is the conditional REML estimate of
$\sigma^2$ at $\widehat{\mathbf{\theta}}_R$, $$\label{eq:REMLsigmasq}
  \widehat{\sigma^2_R}=r^2_{\widehat\theta_R}/(n-p) .$$ It is not
entirely clear how one would define a "REML estimate" of $\mathbf{\beta}$
because the REML criterion, $d_R(\mathbf{\theta},\sigma|\vec
y)$, defined in ([\[eq:REMLdev\]](#eq:REMLdev){reference-type="ref"
reference="eq:REMLdev"}), does not depend on $\mathbf{\beta}$. However, it is
customary (and not unreasonable) to use
$\widehat{\mathbf{\beta}}_R=\widehat{\mathbf{\beta}}_{\widehat{\mathbf{\theta}}_R}$ as
the REML estimate of $\mathbf{\beta}$.

Step-by-step Evaluation of the Profiled Deviance {#sec:stepByStep}
------------------------------------------------

An object returned by is an S4-classed object [@R:Chambers:2008] of
class . A special utility function, , creates a function to evaluate the
deviance from such an object. Because the deviance evaluation function
may be called many, many times for different values of the model
parameters and may need to handle large data structures very
efficiently, the default deviance evaluator uses compiled code, that is
written in C++ using the facilities of the package [@Rcpp]. As such,
this function is not very illuminating

unless you are willing to go digging into the source code for the
compiled function called .

However, if we set the optional argument and suppress the compiled
deviance evaluation then the deviance function follows the deviance
evaluation in more easily understood constructions.

Providing for an -based evaluation of the deviance, in addition to the
compiled code, allows us to check for consistent results. For example,
suppose that we wished to check that the compiled and -based deviance
functions agreed at the initial values of the parameter
$\mathbf{\theta}=[1,0,1]'$

Without going into detail let us just point out that, because the
function is created within another function, , has access to the
structures in the model, , which internally is identified as . The
deviance evaluation uses the value of $\mathbf{\theta}$, called , and
information in three of the model slots: the random-effects structure,
called , the fixed-effects structure, called , and the response
structure, called . There is also a diagonal matrix of weights, called ,
available to this function but for our evaluation that matrix is the
identity and we will ignore it.

The function, with line numbers, is

Lines 3 to 9 create local versions of $\Lambda$ and $\mathbf{L}$ from the
information in the slot and the argument .

In lines 3 to 5 the argument is checked for correct mode and length and
whether it violates the lower bounds stored in the slot.

In lines 6 and 7 the current value of $\Lambda$ is created from a
template, stored in the slot, the value of and the index vector, which
maps elements of $\mathbf{\theta}$ to the non-zero elements of $\Lambda$.

As mentioned above, is the identity in this case so line 8 amounts to
assigning the matrix in the slot to the local variable . Then line 9
evaluates the sparse Cholesky factor, , by updating the template in the
slot. The optional argument, , to the method is the multiple of the
identity matrix to add to the of the second argument. This produces the
factor $\mathbf{L}_\theta$ defined in
([\[eq:sparseCholeskyP\]](#eq:sparseCholeskyP){reference-type="ref"
reference="eq:sparseCholeskyP"}).

Lines 10 to 14 produce $\mathbf{Z}'\mathbf{y}$, $\mathbf{X}'\mathbf{y}$ and
$\mathbf{Z}'\mathbf{X}$, taking into account possible weights, from the
(square root of the residual weights) slot, and/or an offset, from the
slot. To avoid confusion when incorporating weights, these values are
stored as , and .

Because the -based evaluation function, , must be able to handle
non-default values of arguments such as and to the function, we will
skip over some of the details and concentrate on the parts that
correspond to formulas in the previous section.

First we derive the matrix $\Lambda_\theta$ from $\mathbf{\theta}$. A
template version of $\Lambda$ is stored in the random-effects structure
slot, called , of the model object. The slot contains an index vector
indicating which element of is to be used to replace an element of the
slot in

object [@R:Chambers:2008] containing a set of well-defined "slots". For
evaluation of the devai environment, accessed with the extractor. This
environment contains several matrices and vectors that are used in the
evaluation of the profiled deviance. In this section we use these
matrices and vectors from one of our examples to explicitly trace the
steps in evaluating the profiled deviance. This level of detail is
provided for those whose style of learning is more of a "hands on" style
and for those who may want to program modifications of this approach.

Consider our model , fit as

The environment of the model contains the converged parameter vector,
$\mathbf{\theta}$ (), the relative covariance factor, $\Lambda_\theta$ (),
the sparse Cholesky factor, $\vec
L_\theta$ (), the matrices $\mathbf{R}_{ZX}$ () and $\mathbf{R}_X$ (), the
conditional mode, $\tilde{\mathbf{u}}$ (), and the conditional estimate,
$\widehat{\mathbf{\beta}}_\theta$ (). The permutation represented by $\mathbf{P}$
is contained in the sparse Cholesky representation, .

Although the model matrices, $\mathbf{X}$ () and $\mathbf{Z}'$ (), and the
response vector, $\mathbf{y}_{\text{obs}}$ (), are available in the
environment, many of the products that involve only these fixed values
are precomputed and stored separately under the names
($\mathbf{X}'\mathbf{X}$), , and .

To provide easy access to the objects in the environment of we attach it
to the search path.

Please note that this is done here for illustration only. The practice
of attaching a list or a data frame or, less commonly, an environment in
an session is overused, somewhat dangerous (because of the potential of
forgetting to detach it later) and discouraged. The preferred practice
is to use the function to gain access by name to components of such
composite objects. For this section of code, however, using or would
quickly become very tedious and we use instead.

To update the matrix $\Lambda_\theta$ to a new value of $\mathbf{\theta}$ we
need to know which of the non-zeros in $\Lambda$ are updated from which
elements of $\mathbf{\theta}$. Recall that the dimension of $\mathbf{\theta}$ is
small (3, in this case) but $\Lambda$ is potentially large ($18\times18$
with $54$ non-zeros). The environment contains an integer vector that
maps the elements of to the non-zeros in .

Suppose we wish to recreate the evaluation of the profiled deviance at
the initial value of $\mathbf{\theta}=(1,0,1)$. We begin by updating
$\Lambda_\theta$ and forming the product
$\mathbf{u}'=\Lambda_\theta'\mathbf{Z}'$

The Cholesky factor object, , can be updated from without forming
$\mathbf{u}'\mathbf{u}+\mathbf{I}$ explicitly. The optional argument to the
method specifies a multiple of the identity to be added to
$\mathbf{u}'\mathbf{u}$

Then we evaluate and according to
([\[eq:RZXdef\]](#eq:RZXdef){reference-type="ref"
reference="eq:RZXdef"}) and
([\[eq:RXdef\]](#eq:RXdef){reference-type="ref" reference="eq:RXdef"})

Solving ([\[eq:bigPLS\]](#eq:bigPLS){reference-type="ref"
reference="eq:bigPLS"}) for $\tilde{\mathbf{u}}$ and
$\widehat{\mathbf{\beta}}_\theta$ is done in stages. Writing $\mathbf{c}_u$ and
$\mathbf{c}_\beta$ for the intermediate results that satisfy
$$\label{eq:stage1}
    \begin{bmatrix}
    \mathbf{L}_\theta & \mathbf{0}\\
    \mathbf{R}_{ZX}' & \mathbf{R}_X'
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{c}_u\\\mathbf{c}_\beta
  \end{bmatrix}=
  \begin{bmatrix}
    \mathbf{P}\Lambda_\theta'\mathbf{Z}'\mathbf{y}_{\text{obs}}\\
    \mathbf{X}'\mathbf{y}_{\text{obs}} .
  \end{bmatrix}$$ we evaluate

The next set of equations to solve is $$\label{eq:stage2}
    \begin{bmatrix}
    \mathbf{L}_\theta' & \mathbf{R}_{ZX}\\
    \mathbf{0} & \mathbf{R}_X
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{P}\tilde{\mathbf{u}}\\
    \widehat{\mathbf{\beta}}_\theta
  \end{bmatrix}=
  \begin{bmatrix}
    \mathbf{c}_U\\\mathbf{c}_\beta
  \end{bmatrix}.$$

We can now create the conditional mean, , the penalized residual sum of
squares, , the logarithm of the square of the determinant of $\mathbf{L}$, ,
and the profiled deviance, which, fortuitously, equals the value shown
earlier.

to avoid later name clashes.

In terms of the calculations performed, these steps describe exactly the
evaluation of the profiled deviance in . The actual function for
evaluating the deviance, accessible as , ---FIXME--- is a slightly
modified version of what is shown above. However, the modifications are
only to avoid creating copies of potentially large objects and to allow
for cases where the model matrix, $\mathbf{X}$, is sparse. In practice,
unless the optional argument is given, the profiled deviance is
evaluated in compiled code, providing a speed boost, but the code can be
used if desired. This allows for checking the results from the compiled
code and can also be used as a template for extending the computational
methods to other types of models.

Generalizing to Other Forms of Mixed Models {#sec:general}
-------------------------------------------

In later chapters we cover the theory and practice of generalized linear
mixed models (GLMMs), nonlinear mixed models (NLMMs) and generalized
nonlinear mixed models (GNLMMs). Because quite a bit of the theoretical
and computational methodology covered in this chapter extends to those
models we will cover the common aspects here.

### Descriptions of the Model Forms {#sec:modelForms}

We apply the name "generalized" to models in which the conditional
distribution, $(\mathcal{Y}|\mathcal{U}=\mathbf{u})$, is not required to be Gaussian but
does preserve some of the properties of the spherical Gaussian
conditional distribution $$(\mathcal{Y}|\mathcal{U}=\mathbf{u})\sim\mathcal{N}
  (\mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\mathbf{\beta},\sigma^2\mathbf{I}_n)$$ from
the linear mixed model. In particular, the components of $\mathcal{Y}$ are
*conditionally independent*, given $\mathcal{U}=\mathbf{u}$. Furthermore, $\mathbf{u}$
affects the distribution only through the conditional mean, which we
will continue to write as $\mathbf{\mu}$, and it affects the conditional mean
only through the linear predictor,
$\mathbf{\gamma}=\mathbf{Z}\Lambda_\theta\mathbf{u}+\mathbf{X}\mathbf{\beta}$.

Typically we do not have $\mathbf{\mu}=\mathbf{\gamma}$, however. The elements of
the linear predictor, $\mathbf{\gamma}$, can be positive or negative or zero.
Theoretically they can take on any value between $-\infty$ and $\infty$.
But many distributional forms used in GLMMs put constraints on the value
of the mean. For example, the mean of a Bernoulli random variable,
modeling a binary response, must be in the range $0<\mu<1$ and the mean
of a Poisson random variable, modeling a count, must be positive. To
achieve these constraints we write the conditional mean, $\mathbf{\mu}$, as a
transformation of the unbounded predictor, written $\vec\eta$. For
historical, and some theoretical, reasons the inverse of this
transformation is called the *link function*, written
$$\label{eq:linkfunction}
  \vec\eta=\mathbf{g}(\mathbf{\mu}) ,$$ and the transformation we want is called
the *inverse link*, written $\mathbf{g}^{-1}$.

Both $\mathbf{g}$ and $\mathbf{g}^{-1}$ are determined by scalar functions, $g$
and $g^{-1}$, respectively, applied to the individual components of the
vector argument. That is, $\vec\eta$ must be $n$-dimensional and the
vector-valued function $\mathbf{\mu}=\mathbf{g}^{-1}(\vec\eta)$ is defined by the
component functions $\mu_i=g^{-1}(\eta_i),\,i=1,\dots,n$. Among other
things, this means that the Jacobian matrix of the inverse link,
$\frac{d\mathbf{\mu}}{d\vec\eta}$, will be diagonal.

Because the link function, $\mathbf{g}$, and the inverse link, $\vec
g^{-1}$, are nonlinear functions (there would be no purpose in using a
linear link function) many people use the terms "generalized linear
mixed model" and "nonlinear mixed model" interchangeably. We reserve the
term "nonlinear mixed model" for the type of models used, for example,
in pharmacokinetics and pharmacodynamics, where the conditional
distribution is a spherical multivariate Gaussian
$$\label{eq:NLMMconddist}
  (\mathcal{Y}|\mathcal{U}=\mathbf{u})\sim\mathcal{N}(\mathbf{\mu}, \sigma^2\mathbf{I}_n)$$ but
$\mathbf{\mu}$ depends nonlinearly on $\mathbf{\gamma}$. For NLMMs the length of
the linear predictor, $\mathbf{\gamma}$, is a multiple, $ns$, of $n$, the
length of $\mathbf{\mu}$.

Like the map from $\vec\eta$ to $\mathbf{\mu}$, the map from $\mathbf{\gamma}$ to
$\mathbf{\mu}$ has a "diagonal" property, which we now describe. If we use
$\mathbf{\gamma}$ to fill the columns of an $n\times
s$ matrix, $\Gamma$, then $\mu_i$ depends only on the $i$th row of
$\Gamma$. In fact, $\mu_i$ is determined by a nonlinear model function,
$f$, applied to the $i$ row of $\Gamma$. Writing
$\mathbf{\mu}=\mathbf{f}(\mathbf{\gamma})$ based on the component function $f$, we see
that the Jacobian of $\mathbf{f}$, $\frac{d\mathbf{\mu}}{d\mathbf{\gamma}}$, will be
the vertical concatenation of $s$ diagonal $n\times n$ matrices.

Because we will allow for generalized nonlinear mixed models (GNLMMs),
in which the mapping from $\mathbf{\gamma}$ to $\mathbf{\mu}$ has the form
$$\label{eq:GammaEtaMu}
  \mathbf{\gamma}\;\rightarrow\;\vec\eta\;\rightarrow\;\mathbf{\mu} ,$$ we will
use ([\[eq:GammaEtaMu\]](#eq:GammaEtaMu){reference-type="ref"
reference="eq:GammaEtaMu"}) in our definitions.

### Determining the Conditional Mode, $\tilde{\mathbf{u}}$ {#sec:conditionalmode}

For all these types of mixed models, the conditional distribution,
$(\mathcal{U}|\mathcal{Y}=\mathbf{y}_{\text{obs}})$, is a continuous distribution for
which we can determine the unscaled conditional density, $h(\mathbf{u})$. As
for linear mixed models, we define the conditional mode,
$\tilde{\mathbf{u}}$, as the value that maximizes the unscaled conditional
density.

Determining the conditional mode, $\tilde{\mathbf{u}}$, in a nonlinear mixed
model is a penalized nonlinear least squares (PNLS) problem
$$\label{eq:PNLSprob}
  \tilde{\mathbf{u}}=\arg\min_{\mathbf{u}}\|\mathbf{y}_{\text{obs}}-\mathbf{\mu}\|^2+\|\mathbf{u}\|^2$$
which we solve by adapting the iterative techniques, such as the
Gauss-Newton method [@bateswatts88:_nonlin Sect. 2.2.1], used for
nonlinear least squares. Starting at an initial value, $\vec
u^{(0)}$, (the bracketed superscript denotes the iteration number) with
conditional mean, $\mathbf{\mu}^{(0)}$, we determine an increment
$\vec\delta^{(1)}$ by solving the penalized linear least squares
problem, $$\label{eq:PNLSiter}
  \vec\delta^{(1)}=\arg\min_{\vec\delta}\left\|
    \begin{bmatrix}
      \mathbf{y}_{\text{obs}}-\mathbf{\mu}^{(0)}\\
      \mathbf{0}-\mathbf{u}^{(0)}
    \end{bmatrix} -
    \begin{bmatrix}
      \mathbf{u}^{(0)}\\
      \mathbf{I}_q
    \end{bmatrix}\vec\delta\right\|^2$$ where $$\label{eq:Udef}
  \mathbf{u}^{(0)}=\left.\frac{d\mathbf{\mu}}{d\mathbf{u}}\right|_{\mathbf{u}^{(0)}} .$$
Naturally, we use the sparse Cholesky decomposition, $\vec
L_\theta^{(0)}$, satisfying $$\label{eq:sparseCholPNLS}
  \mathbf{L}_\theta^{(0)}\left(\mathbf{L}_\theta^{(0)}\right)=
  \mathbf{P}\left[\left(\mathbf{u}^{(0)}\right)'\mathbf{u}^{(0)}+\mathbf{I}_q\right]
  \mathbf{P}'$$ to determine this increment. The next iteration begins
at $$\label{eq:u1}
  \mathbf{u}^{(1)} = \mathbf{u}^{(0)} + k \vec\delta^{(1)}$$ where $k$ is the
step factor chosen, perhaps by step-halving [@bateswatts88:_nonlin
Sect. 2.2.1], to ensure that the penalized residual sum of squares
decreases at each iteration. Convergence is declared when the
orthogonality convergence criterion [@bateswatts88:_nonlin Sect. 2.2.3]
is below some pre-specified tolerance.

The *Laplace approximation* to the deviance is $$\label{eq:Laplace}
  d(\mathbf{\theta},\mathbf{\beta},\sigma|\mathbf{y}_{\text{obs}})\approx
  n\log(2\pi\sigma^2)+2\log|\mathbf{L}_{\theta,\beta}|+
  \frac{r^2_{\theta,\beta}}{\sigma^2},$$ where the Cholesky factor,
$\mathbf{L}_{\theta,\beta}$, and the penalized residual sum of squares,
$r^2_{\theta,\beta}$, are both evaluated at the conditional mode,
$\tilde{\mathbf{u}}$. The Cholesky factor depends on $\mathbf{\theta}$,
$\mathbf{\beta}$ and $\mathbf{u}$ for these models but typically the dependence
on $\mathbf{\beta}$ and $\mathbf{u}$ is weak.

Chapter Summary {#sec:lmmsummary}
---------------

The definitions and the computational results for maximum likelihood
estimation of the parameters in linear mixed models were summarized in .
A key computation is evaluation of the sparse Cholesky factor,
$\Lambda_\theta$, satisfying
([\[eq:sparseCholeskyP\]](#eq:sparseCholeskyP){reference-type="ref"
reference="eq:sparseCholeskyP"}), $$\mathbf{L}_\theta\mathbf{L}_\theta'
   = \mathbf{P}\left(\Lambda_\theta'\mathbf{Z}'\mathbf{Z}\Lambda_\theta
   +\mathbf{I}_q\right)\mathbf{P}',$$ where $\mathbf{P}$ represents the
fill-reducing permutation determined during the symbolic phase of the
sparse Cholesky decomposition.

An extended decomposition
([\[eq:bigdecomp\]](#eq:bigdecomp){reference-type="ref"
reference="eq:bigdecomp"}) provides the $q\times p$ matrix $\mathbf{R}_{ZX}$
and the $p\times p$ upper triangular $\mathbf{R}_X$ that are used to
determine the conditional mode $\tilde{\mathbf{u}}_\theta$, the conditional
estimate $\widehat{\mathbf{\beta}}_\theta$, and the minimum penalized
residual sum of squares, $r^2_\theta$, from which the profiled deviance
$$\tilde{d}(\mathbf{\theta}|\mathbf{y}_{\text{obs}})
  =2\log|\mathbf{L}_\theta|+n\left[1 +
    \log\left(\frac{2 \pi r^2_\theta}{n}\right)\right]$$ or the profile
REML criterion $$\tilde{d}_R(\mathbf{\theta}|\mathbf{y}_{\text{obs}})=
  2\log\left(|\mathbf{L}_\theta||\mathbf{R}_X|\right)+(n-p)
  \left[1+\log\left(\frac{2\pi r^2_\theta}{n-p}\right)\right]$$ can be
evaluated and optimized (minimized) with respect to $\mathbf{\theta}$.
