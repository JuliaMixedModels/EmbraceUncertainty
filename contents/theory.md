Computational Methods for Mixed Models {#chap:computational}
======================================

In this chapter we describe some of the details of the computational
methods for fitting linear mixed models, as implemented in the package,
and the theoretical development behind these methods. We also provide
the basis for later generalizations to models for non-Gaussian responses
and to models in which the relationship between the conditional mean,
$\vec\mu$, and the linear predictor,
$\vec\gamma=\vec X\vec\beta+\vec Z\vec b=
\vec Z\Lambda_\theta\vec u+\vec X\vec\beta$, is a nonlinear
relationship.

This material is directed at those readers who wish to follow the theory
and methodology of linear mixed models and how both can be extended to
other forms of mixed models. Readers who are less interested in the
"how" and the "why" of fitting mixed models than in the results
themselves should not feel obligated to master these details.

We begin by reviewing the definition of linear mixed-effects models and
some of the basics of the computational methods, as given in .

Definitions and Basic Results {#sec:defnLMM}
-----------------------------

As described in , a linear mixed-effects model is based on two
vector-valued random variables: the $q$-dimensional vector of random
effects, $\bc B$, and the $n$-dimensional response vector, $\bc Y$.
Equation ([\[eq:LMMdist\]](#eq:LMMdist){reference-type="ref"
reference="eq:LMMdist"}) defines the unconditional distribution of
$\bc B$ and the conditional distribution of $\bc Y$, given
$\bc B=\vec b$, as multivariate Gaussian distributions of the form
$$\begin{aligned}
    (\bc Y|\bc B=\vec b)&\sim\mathcal{N}(\vec X\vec\beta+\vec Z\vec
    b,\sigma^2\vec I)\\
    \bc{B}&\sim\mathcal{N}(\vec0,\Sigma_\theta) .
  \end{aligned}$$

The $q\times q$, symmetric, variance-covariance matrix,
$\mathrm{Var}(\bc B)=\Sigma_\theta$, depends on the *variance-component
parameter vector*, $\vec\theta$, and is *positive semidefinite*, which
means that $$\label{eq:posSemiDef}
  \vec b\trans\Sigma_\theta\vec b\ge0,\quad\forall\,\vec b\ne\vec 0 .$$
(The symbol $\forall$ denotes "for all".) The fact that $\Sigma_\theta$
is positive semidefinite does not guarantee that $\Sigma_\theta^{-1}$
exists. We would need a stronger property, $\vec
b\trans\Sigma_\theta\vec b>0,\,\forall\,\vec b\ne\vec 0$, called
positive definiteness, to ensure that $\Sigma_\theta^{-1}$ exists.

Many computational formulas for linear mixed models are written in terms
of $\Sigma_\theta^{-1}$. Such formulas will become unstable as
$\Sigma_\theta$ approaches singularity. And it can do so. It is a fact
that singular (i.e. non-invertible) $\Sigma_\theta$ can and do occur in
practice, as we have seen in some of the examples in earlier chapters.
Moreover, during the course of the numerical optimization by which the
parameter estimates are determined, it is frequently the case that the
deviance or the REML criterion will need to be evaluated at values of
$\vec\theta$ that produce a singular $\Sigma_\theta$. Because of this we
will take care to use computational methods that can be applied even
when $\Sigma_\theta$ is singular and are stable as $\Sigma_\theta$
approaches singularity.

As defined in ([\[eq:relcovfac\]](#eq:relcovfac){reference-type="ref"
reference="eq:relcovfac"}) a relative covariance factor,
$\Lambda_\theta$, is any matrix that satisfies
$$\Sigma_\theta=\sigma^2\Lambda_\theta\Lambda_\theta\trans .$$ According
to this definition, $\Sigma$ depends on both $\sigma$ and $\theta$ and
we should write it as $\Sigma_{\sigma,\theta}$. However, we will blur
that distinction and continue to write $\text{Var}(\bc
B)=\Sigma_\theta$. Another technicality is that the *common scale
parameter*, $\sigma$, can, in theory, be zero. We will show that in
practice the only way for its estimate, $\widehat{\sigma}$, to be zero
is for the fitted values from the fixed-effects only, $\vec
X\widehat{\vec\beta}$, to be exactly equal to the observed data. This
occurs only with data that have been (incorrectly) simulated without
error. In practice we can safely assume that $\sigma>0$. However,
$\Lambda_\theta$, like $\Sigma_\theta$, can be singular.

Our computational methods are based on $\Lambda_\theta$ and do not
require evaluation of $\Sigma_\theta$. In fact, $\Sigma_\theta$ is
explicitly evaluated only at the converged parameter estimates.

The spherical random effects, $\bc U\sim\mathcal{N}(\vec
0,\sigma^2\vec I_q)$, determine $\bc B$ as $$\label{eq:sphericalRE}
  \bc B=\Lambda_\theta\bc U .$$ Although it may seem more intuitive to
write $\bc U$ as a linear transformation of $\bc B$, we cannot do that
when $\Lambda_\theta$ is singular, which is why
([\[eq:sphericalRE\]](#eq:sphericalRE){reference-type="ref"
reference="eq:sphericalRE"}) is in the form shown.

We can easily verify that
([\[eq:sphericalRE\]](#eq:sphericalRE){reference-type="ref"
reference="eq:sphericalRE"}) provides the desired distribution for
$\bc B$. As a linear transformation of a multivariate Gaussian random
variable, $\bc B$ will also be multivariate Gaussian. Its mean and
variance-covariance matrix are straightforward to evaluate,
$$\label{eq:EB}
  \mathrm{E}[\bc B] = \Lambda_\theta\mathrm{E}[\bc U]=\Lambda_\theta\vec0=\vec0$$
and $$\mathrm{Var}(\bc B)
  \begin{aligned}[t]
    &=\mathrm{E}\left[(\bc B-\mathrm{E}[\bc B])
      (\bc B-\mathrm{E}[\bc B])\trans\right]
    =\mathrm{E}\left[\bc B\bc B\trans\right]\\
    &=\mathrm{E}\left[\Lambda_\theta\,\bc U\bc U\trans\Lambda_\theta\trans\right]
    =\Lambda_\theta\,\mathrm{E}[\bc U\bc U\trans]\Lambda_\theta\trans
    =\Lambda_\theta\,\mathrm{Var}(\bc U)\Lambda_\theta\trans\\
    &=\Lambda_\theta\,\sigma^2\vec I_q\,\Lambda_\theta\trans
    =\sigma^2\Lambda_\theta\Lambda_\theta\trans
    =\Sigma_\theta
  \end{aligned}$$ and have the desired form.

Just as we concentrate on how $\vec\theta$ determines $\Lambda_\theta$,
not $\Sigma_\theta$, we will concentrate on properties of $\bc U$ rather
than $\bc B$. In particular, we now define the model according to the
distributions $$\label{eq:condYgivenU}
  \begin{aligned}
  (\bc Y|\bc U=\vec u)&\sim\mathcal{N}(\vec Z\Lambda_\theta\vec
  u+\vec X\beta,\sigma^2\vec I_n)\\
  \bc U&\sim\mathcal{N}(\vec0,\sigma^2\vec I_q) .
  \end{aligned}$$

To allow for extensions to other types of mixed models we distinguish
between the *linear predictor* $$\label{eq:linearpred}
  \vec\gamma = \vec Z\Lambda_\theta\vec u+\vec X\beta$$ and the
*conditional mean* of $\bc Y$, given $\bc U=\vec u$, which is
$$\label{eq:conditionalMean}
  \vec\mu = \mathrm{E}\left[\bc Y|\bc U=\vec u\right] .$$ For a linear
mixed model $\vec\mu=\vec\gamma$. In other forms of mixed models the
conditional mean, $\vec\mu$, can be a nonlinear function of the linear
predictor, $\vec\gamma$. For some models the dimension of $\vec\gamma$
is a multiple of $n$, the dimension of $\vec\mu$ and $\vec y$, but for a
linear mixed model the dimension of $\vec\gamma$ must be $n$. Hence, the
model matrix $\vec Z$ must be $n\times q$ and $\vec X$ must be
$n\times p$.

The Conditional Distribution $(\mathcal{U}|\mathcal{Y}=\vec y)$ {#sec:conddistUgivenY}
---------------------------------------------------------------

In this chapter it will help to be able to distinguish between the
observed response vector and an arbitrary value of $\bc Y$. For this
chapter only we will write the observed data vector as $\vec
y_{\text{obs}}$, with the understanding that $\vec y$ without the
subscript will refer to an arbitrary value of the random variable $\bc
Y$.

The likelihood of the parameters, $\vec\theta$, $\vec\beta$, and
$\sigma$, given the observed data, $\vec y_{\text{obs}}$, is the
probability density of $\bc Y$, evaluated at $\vec y_{\text{obs}}$.
Although the numerical values of the probability density and the
likelihood are identical, the interpretations of these functions are
different. In the density we consider the parameters to be fixed and the
value of $\vec y$ as varying. In the likelihood we consider $\vec y$ to
be fixed at $\vec y_{\text{obs}}$ and the parameters, $\vec\theta$,
$\vec\beta$ and $\sigma$, as varying.

The natural approach for evaluating the likelihood is to determine the
marginal distribution of $\bc Y$, which in this case amounts to
determining the marginal density of $\bc Y$, and evaluate that density
at $\vec
y_{\text{obs}}$. To follow this course we would first determine the
joint density of $\bc U$ and $\bc Y$, written
$f_{\bc U,\bc Y}(\vec u,\vec y)$, then integrate this density with
respect to $\vec u$ to create the marginal density, $f_{\bc Y}(\vec y)$,
and finally evaluate this marginal density at $\vec y_{\text{obs}}$.

To allow for later generalizations we will change the order of these
steps slightly. We evaluate the joint density function, $f_{\bc U,\bc
  Y}(\vec u,\vec y)$, at $\vec y_{\text{obs}}$, producing the
*unnormalized conditional density*, $h(\vec u)$. We say that $h$ is
"unnormalized" because the conditional density is a multiple of $h$
$$\label{eq:conddenUgivenY}
  f_{\bc U|\bc Y}(\vec u|\vec y_{\text{obs}})=\frac
  {h(\vec u)}{\int_{\mathbb{R}^q}h(\vec u)\,d\vec u}  .$$ In some
theoretical developments the normalizing constant, which is the integral
in the denominator of an expression like
([\[eq:conddenUgivenY\]](#eq:conddenUgivenY){reference-type="ref"
reference="eq:conddenUgivenY"}), is not of interest. Here it is of
interest because the normalizing constant is exactly the likelihood that
we wish to evaluate, $$\label{eq:LMMlikelihood}
  L(\vec\theta,\vec\beta,\sigma|\vec y_{\text{obs}}) =
  \int_{\mathbb{R}^q}h(\vec u)\,d\vec u .$$

For a linear mixed model, where all the distributions of interest are
multivariate Gaussian and the conditional mean, $\vec\mu$, is a linear
function of both $\vec u$ and $\vec\beta$, the distinction between
evaluating the joint density at $\vec y_{\text{obs}}$ to produce
$h(\vec u)$ then integrating with respect to $\vec u$, as opposed to
first integrating the joint density then evaluating at $\vec
y_{\text{obs}}$, is not terribly important. For other mixed models this
distinction can be important. In particular, generalized linear mixed
models, described in
Chap. [\[chap:GLMMbinomial\]](#chap:GLMMbinomial){reference-type="ref"
reference="chap:GLMMbinomial"}, are often used to model a discrete
response, such as a binary response or a count, leading to a joint
distribution for $\bc Y$ and $\bc U$ that is discrete with respect to
one variable, $\vec y$, and continuous with respect to the other,
$\vec u$. In such cases there isn't a joint density for $\bc Y$ and
$\bc U$. The necessary distribution theory for general $\vec y$ and
$\vec u$ is well-defined but somewhat awkward to describe. It is much
easier to realize that we are only interested in the observed response
vector, $\vec y_{\text{obs}}$, not some arbitrary value of $\vec y$, so
we can concentrate on the conditional distribution of $\bc U$ given
$\bc Y=\vec y_{\text{obs}}$. For all the mixed models we will consider,
the conditional distribution, $(\bc
U|\bc Y=\vec y_{\text{obs}})$, is continuous and both the conditional
density, $f_{\bc U|\bc Y}(\vec u|\vec y_{\text{obs}})$, and its
unnormalized form, $h(\vec u)$, are well-defined.

Integrating $h(\vec u)$ in the Linear Mixed Model {#sec:IntegratingH}
-------------------------------------------------

The integral defining the likelihood in
([\[eq:LMMlikelihood\]](#eq:LMMlikelihood){reference-type="ref"
reference="eq:LMMlikelihood"}) has a closed form in the case of a linear
mixed model but not for some of the more general forms of mixed models.
To motivate methods for approximating the likelihood in more general
situations, we describe in some detail how the integral can be evaluated
using the sparse Cholesky factor, $\vec L_\theta$, and the conditional
mode, $$\label{eq:condMode}
  \tilde{\vec u}=\arg\max_{\vec u} f_{\bc U|\bc Y}(\vec u|\vec y_{\text{obs}})=
  \arg\max_{\vec u} h(\vec u) = \arg\max_{\vec u}
  f_{\bc Y|\bc U}(\vec y_{\text{obs}}|\vec u)\,f_{\bc U}(\vec u).$$ The
notation $\arg\max_{\vec u}$ means that $\tilde{\vec u}$ is the value of
$\vec u$ that maximizes the expression that follows.

In general, the *mode* of a continuous distribution is the value of the
random variable that maximizes the density. The value $\tilde{\vec u}$
is called the conditional mode of $\vec u$, given
$\bc Y=\vec y_{\text{obs}}$, because $\tilde{\vec u}$ maximizes the
conditional density of $\bc U$ given $\bc Y=\vec y_{\text{obs}}$. The
location of the maximum can be determined by maximizing the unnormalized
conditional density because $h(\vec u)$ is just a constant multiple of
$f_{\bc U|\bc Y}(\vec u|\vec y_{\text{obs}})$. The last part of
([\[eq:condMode\]](#eq:condMode){reference-type="ref"
reference="eq:condMode"}) is simply a re-expression of $h(\vec u)$ as
the product of $f_{\bc Y|\bc U}(\vec y_{\text{obs}}|\vec u)$ and
$f_{\bc U}(\vec u)$. For a linear mixed model these densities are
$$\begin{aligned}
  \label{eq:densYgivenUandU}
  f_{\bc Y|\bc U}(\vec y|\vec u)&=
  \frac{1}{\left(2\pi\sigma^2\right)^{n/2}}
  \exp\left(-\frac{\left\|\vec y-\vec X\vec\beta-\vec Z\Lambda_\theta\vec u\right\|^2}{2\sigma^2}\right)\\
  f_{\bc U}(\vec u)&=
  \frac{1}{\left(2\pi\sigma^2\right)^{q/2}}\exp\left(-\frac{\|\vec u\|^2}
    {2\sigma^2}\right)\end{aligned}$$ with product $$\label{eq:hudef}
  h(\vec u)=\frac{1}{\left(2\pi\sigma^2\right)^{(n+q)/2}}
  \exp\left(-\frac{\left\|\vec y_{\text{obs}}-\vec X\vec\beta-\vec Z\Lambda_\theta\vec u\right\|^2+\|\vec u\|^2}{2\sigma^2}\right) .$$
On the deviance scale we have $$\label{eq:devh}
  -2\log\left(h(\vec u)\right)=(n+q)\log(2\pi\sigma^2)
  +\frac{\left\|\vec y_{\text{obs}}-\vec X\vec\beta-\vec Z\Lambda_\theta\vec
      u\right\|^2+\|\vec u\|^2}{\sigma^2} .$$ Because
([\[eq:devh\]](#eq:devh){reference-type="ref" reference="eq:devh"})
describes the negative log density, $\tilde{\vec u}$ will be the value
of $\vec u$ that minimizes the expression on the right hand side of
([\[eq:devh\]](#eq:devh){reference-type="ref" reference="eq:devh"}).

The only part of the right hand side of
([\[eq:devh\]](#eq:devh){reference-type="ref" reference="eq:devh"}) that
depends on $\vec u$ is the numerator of the second term. Thus
$$\label{eq:PLSsol}
  \tilde{\vec u}=\arg\min_{\vec u} \left\|
    \vec y_{\text{obs}}-\vec X\vec\beta-\vec Z\Lambda_\theta\vec u\right\|^2+
  \|\vec u\|^2.$$ The expression to be minimized, called the *objective
function*, is described as a *penalized residual sum of squares* (PRSS)
and the minimizer, $\tilde{\vec u}$, is called the *penalized least
squares* (PLS) solution. They are given these names because the first
term in the objective, $\left\| \vec y_{\text{obs}}-\vec X\vec\beta-\vec
  Z\Lambda_\theta\vec u\right\|^2$, is a sum of squared residuals, and
the second term, $\|\vec u\|^2$, is a penalty on the length, $\|\vec
u\|$, of $\vec u$. Larger values of $\vec u$ (in the sense of greater
lengths as vectors) incur a higher penalty.

The PRSS criterion determining the conditional mode balances fidelity to
the observed data (i.e. producing a small residual sum of squares)
against simplicity of the model (small $\|\vec u\|$). We refer to this
type of criterion as a smoothing objective, in the sense that it seeks
to smooth out the fitted response by reducing model complexity while
still retaining reasonable fidelity to the observed data.

For the purpose of evaluating the likelihood we will regard the PRSS
criterion as a function of the parameters, given the data, and write its
minimum value as $$\label{eq:r2thetabeta}
  r^2_{\theta,\beta}=\min_{\vec u} \left\|
    \vec y_{\text{obs}}-\vec X\vec\beta-\vec Z\Lambda_\theta\vec u\right\|^2+ \|\vec u\|^2.$$
Notice that $\vec\beta$ only enters the right hand side of
([\[eq:r2thetabeta\]](#eq:r2thetabeta){reference-type="ref"
reference="eq:r2thetabeta"}) through the linear predictor expression. We
will see that $\tilde{\vec u}$ can be determined by a direct (i.e.
non-iterative) calculation and, in fact, we can minimize the PRSS
criterion with respect to $\vec u$ and $\vec\beta$ simultaneously
without iterating. We write this minimum value as $$\label{eq:r2theta}
  r^2_\theta=\min_{\vec u,\vec\beta} \left\|
    \vec y_{\text{obs}}-\vec X\vec\beta-\vec Z\Lambda_\theta\vec u\right\|^2+ \|\vec u\|^2.$$
The value of $\vec\beta$ at the minimum is called the conditional
estimate of $\vec\beta$ given $\vec\theta$, written
$\widehat{\vec\beta}_\theta$.

Determining the PLS Solutions, $\tilde{\vec u}$ and $\widehat{\vec\beta}_\theta$ {#sec:PLSsol}
--------------------------------------------------------------------------------

One way of expressing a penalized least squares problem like
([\[eq:r2thetabeta\]](#eq:r2thetabeta){reference-type="ref"
reference="eq:r2thetabeta"}) is by incorporating the penalty as
"pseudo-data" in an ordinary least squares problem. We extend the
"response vector", which is $\vec y_{\text{obs}}-\vec X\vec\beta$ when
we minimize with respect to $\vec u$ only, with $q$ responses that are 0
and we extend the predictor expression, $\vec Z\Lambda_\theta\vec u$
with $\vec I_q\vec u$. Writing this as a least squares problem produces
$$\label{eq:PLSLMM}
  \tilde{\vec u}=\arg\min_{\vec u}\left\|
    \begin{bmatrix}
      \vec y_{\text{obs}}-\vec X\vec\beta\\
      \vec 0
    \end{bmatrix} -
    \begin{bmatrix}
      \vec Z\Lambda_\theta\\
      \vec I_q
    \end{bmatrix}\vec u\right\|^2$$ with a solution that satisfies
$$\label{eq:LMMPLSsol}
  \left(\Lambda_\theta\trans\vec Z\trans\vec Z\Lambda_\theta+\vec
    I_q\right)\tilde{\vec u}
  =\Lambda_\theta\trans\vec Z\trans\left(\vec y_{\text{obs}}-\vec X\vec\beta\right).$$

To evaluate $\tilde{\vec u}$ we form the *sparse Cholesky factor*,
$\vec L_\theta$, which is a lower triangular $q\times q$ matrix that
satisfies $$\label{eq:sparseCholesky}
  \vec L_\theta\vec L_\theta\trans=
  \Lambda_\theta\trans\vec Z\trans\vec Z\Lambda_\theta+\vec I_q .$$ The
actual evaluation of the sparse Cholesky factor, $\vec L_\theta$, often
incorporates a *fill-reducing permutation*, which we describe next.

### The Fill-reducing Permutation, $\vec P$ {#sec:fill-reducingP}

In earlier chapters we have seen that often the random effects vector is
re-ordered before $\vec L_\theta$ is created. The re-ordering or
permutation of the elements of $\vec u$ and, correspondingly, the
columns of the model matrix, $\vec Z\Lambda_\theta$, does not affect the
theory of linear mixed models but can have a profound effect on the time
and storage required to evaluate $\vec L_\theta$ in large problems. We
write the effect of the permutation as multiplication by a $q\times q$
*permutation matrix*, $\vec P$, although in practice we apply the
permutation without ever constructing $\vec
P$. That is, the matrix $\vec P$ is a notational convenience only.

The matrix $\vec P$ consists of permuted columns of the identity matrix,
$\vec I_q$, and it is easy to establish that the inverse permutation
corresponds to multiplication by $\vec P\trans$. Because multiplication
by $\vec P$ or by $\vec P\trans$ simply re-orders the components of a
vector, the length of the vector is unchanged. Thus,
$$\label{eq:orthogonalP}
  \|\vec P\vec u\|^2= \|\vec u\|^2 = \|\vec P\trans\vec u\|^2$$ and we
can express the penalty in
([\[eq:r2theta\]](#eq:r2theta){reference-type="ref"
reference="eq:r2theta"}) in any of these three forms. The properties of
$\vec P$ that it preserves lengths of vectors and that its transpose is
its inverse are summarized by stating that $\vec P$ is an *orthogonal
matrix*.

The permutation represented by $\vec P$ is determined from the structure
of $\Lambda_\theta\trans\vec Z\trans\vec
Z\Lambda_\theta+\vec I_q$ for some initial value of $\vec\theta$. The
particular value of $\vec\theta$ does not affect the result because the
permutation depends only the positions of the non-zeros, not the
numerical values at these positions.

Taking into account the permutation, the sparse Cholesky factor, $\vec
L_\theta$, is defined to be the sparse, lower triangular, $q\times q$
matrix with positive diagonal elements satisfying
$$\label{eq:sparseCholeskyP}
  \vec L_\theta\vec L_\theta\trans
   = \vec P\left(\Lambda_\theta\trans\vec Z\trans\vec Z\Lambda_\theta
   +\vec I_q\right)\vec P\trans.$$
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

Many sparse matrix methods, including the sparse Cholesky decomposition,
are performed in two stages: the *symbolic phase* in which the locations
of the non-zeros in the result are determined and the *numeric phase* in
which the numeric values at these positions are evaluated. The symbolic
phase for the decomposition
([\[eq:sparseCholeskyP\]](#eq:sparseCholeskyP){reference-type="ref"
reference="eq:sparseCholeskyP"}), which includes determining the
permutation, $\vec P$, need only be done once. Evaluation of $\vec
L_\theta$ for subsequent values of $\vec\theta$ requires only the
numeric phase, which typically is much faster than the symbolic phase.

The permutation, $\vec P$, serves two purposes. The first, and more
important purpose, is to reduce the number of non-zeros in the factor,
$\vec L_\theta$. The factor is potentially non-zero at every non-zero
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
  \tilde{\vec u}=\arg\min_{\vec u}\left\|
    \begin{bmatrix}
      \vec y_{\text{obs}}-\vec X\vec\beta\\
      \vec 0
    \end{bmatrix} -
    \begin{bmatrix}
      \vec Z\Lambda_\theta\vec P\trans\\
      \vec P\trans
    \end{bmatrix}\vec P\vec u\right\|^2$$ and the system of linear
equations satisfied by $\tilde{\vec u}$ is $$\label{eq:LMMPLSsolP}
  \vec L_\theta\vec L_\theta\trans\vec P\tilde{\vec u}=
  \vec P\left(\Lambda_\theta\trans\vec Z\trans\vec Z\Lambda_\theta+\vec
    I_q\right)\vec P\trans\vec P\tilde{\vec u}
  =\vec P\Lambda_\theta\trans\vec Z\trans\left(\vec y_{\text{obs}}-\vec X\vec\beta\right) .$$

Obtaining the Cholesky factor, $\vec L_\theta$, may not seem to be great
progress toward determining $\tilde{\vec u}$ because we still must solve
([\[eq:LMMPLSsolP\]](#eq:LMMPLSsolP){reference-type="ref"
reference="eq:LMMPLSsolP"}) for $\tilde{\vec u}$. However, it is the key
to the computational methods in the package. The ability to evaluate
$\vec L_\theta$ rapidly for many different values of $\vec\theta$ is
what makes the computational methods in feasible, even when applied to
very large data sets with complex structure. Once we evaluate
$\vec L_\theta$ it is straightforward to solve
([\[eq:LMMPLSsolP\]](#eq:LMMPLSsolP){reference-type="ref"
reference="eq:LMMPLSsolP"}) for $\tilde{\vec u}$ because $\vec L_\theta$
is triangular.

In we will describe the steps in determining this solution but first we
will show that the solution, $\tilde{\vec u}$, and the value of the
objective at the solution, $r^2_{\theta,\beta}$, do allow us to evaluate
the deviance.

### The Value of the Deviance and Profiled Deviance {#sec:solvingtildeu}

After evaluating $\vec L_\theta$ and using that to solve for
$\tilde{\vec u}$, which also produces $r^2_{\beta,\theta}$, we can write
the PRSS for a general $\vec u$ as $$\label{eq:PRSSwithL}
\left\|
    \vec y_{\text{obs}}-\vec X\vec\beta-\vec Z\Lambda_\theta\vec u\right\|^2+ \|\vec u\|^2=
  r^2_{\theta,\beta}+\|\vec L_\theta\trans(\vec u-\tilde{\vec u})\|^2$$
which finally allows us to evaluate the likelihood. We plug the right
hand side of ([\[eq:PRSSwithL\]](#eq:PRSSwithL){reference-type="ref"
reference="eq:PRSSwithL"}) into the definition of $h(\vec u)$ and apply
the change of variable $$\label{eq:changeVar}
  \vec z=\frac{\vec L_\theta\trans(\vec u-\tilde{\vec u})}{\sigma} .$$
The determinant of the Jacobian of this transformation,
$$\label{eq:tranJac}
  \left|\frac{d\vec z}{d\vec u}\right|=
  \left|\frac{\vec L_\theta\trans}{\sigma}\right|=
  \frac{|\vec L_\theta|}{\sigma^q}$$ is required for the change of
variable in the integral. We use the letter $\vec z$ for the transformed
value because we will rearrange the integral to have the form of the
integral of the density of the standard multivariate normal
distribution. That is, we will use the result $$\label{eq:stdnormint}
   \int_{\mathbb{R}^q}\frac{e^{-\|\vec z\|^2/2}}{(2\pi)^{q/2}}\,d\vec
   z = 1.$$

Putting all these pieces together gives $$\label{eq:hintegral}
  \begin{aligned}
    L(\theta,\beta,\sigma)&=\int_{\mathbb{R}^q}h(\vec u)\,d\vec u\\
    &=\int_{\mathbb{R}^q}\frac{1}{(2\pi\sigma^2)^{(n+q)/2}}
    \exp\left(-\frac{r^2_{\theta,\beta}+\|\vec L_\theta\trans(\vec
        u-\tilde{\vec u})\|^2}{2\sigma^2}\right)d\vec u\\
    &=\frac{\exp\left(-\frac{r^2_{\theta,\beta}}{2\sigma^2}\right)}
    {(2\pi\sigma^2)^{n/2}}\int_{\mathbb{R}^q}\frac{1}{(2\pi)^{q/2}}
    \exp\left(-\frac{\|\vec L_\theta\trans(\vec u-\tilde{\vec
          u})\|^2}{2\sigma^2}\right)
    \frac{|\vec L_\theta|}{|\vec L_\theta|}
    \frac{d\vec
      u}{\sigma^q}\\
    &=\frac{\exp\left(-\frac{r^2_{\theta,\beta}}{2\sigma^2}\right)}
    {(2\pi\sigma^2)^{n/2}|\vec L_\theta|}
    \int_{\mathbb{R}^q}\frac{e^{-\|\vec
        z\|^2/2}}{(2\pi)^{q/2}}\,d\vec z\\
    &=\frac{\exp\left(-\frac{r^2_{\theta,\beta}}{2\sigma^2}\right)}
    {(2\pi\sigma^2)^{n/2}|\vec L_\theta|} .
  \end{aligned}$$

The deviance can now be expressed as
$$d(\vec\theta,\vec\beta,\sigma|\vec y_{\text{obs}})=
  -2\log\left(L(\vec\theta,\vec\beta,\sigma|\vec y_{\text{obs}})\right)
  =n\log(2\pi\sigma^2)+2\log|\vec L_\theta|+
  \frac{r^2_{\beta,\theta}}{\sigma^2},$$ as stated in
([\[eq:LMMdeviance\]](#eq:LMMdeviance){reference-type="ref"
reference="eq:LMMdeviance"}). The maximum likelihood estimates of the
parameters are those that minimize this deviance.

Equation ([\[eq:LMMdeviance\]](#eq:LMMdeviance){reference-type="ref"
reference="eq:LMMdeviance"}) is a remarkably compact expression,
considering that the class of models to which it applies is very large
indeed. However, we can do better than this if we notice that
$\vec\beta$ affects
([\[eq:LMMdeviance\]](#eq:LMMdeviance){reference-type="ref"
reference="eq:LMMdeviance"}) only through $r^2_{\beta,\theta}$, and, for
any value of $\vec\theta$, minimizing this expression with respect to
$\vec\beta$ is just an extension of the penalized least squares problem.
Let $\widehat{\vec\beta}_\theta$ be the value of $\vec\beta$ that
minimizes the PRSS simultaneously with respect to $\vec\beta$ and
$\vec u$ and let $r^2_\theta$ be the PRSS at these minimizing values.
If, in addition, we set $\widehat{\sigma^2}_\theta=r^2_\theta/n$, which
is the value of $\sigma^2$ that minimizes the deviance for a given value
of $r^2_\theta$, then the *profiled deviance*, which is a function of
$\vec\theta$ only, becomes $$\label{eq:LMMprofdeviance}
  \tilde{d}(\vec\theta|\vec y_{\text{obs}})
  =2\log|\vec L_\theta|+n\left[1 +
    \log\left(\frac{2 \pi r^2_\theta}{n}\right)\right].$$

Numerical optimization (minimization) of $\tilde{d}(\vec\theta|\vec
y_{\text{obs}})$ with respect to $\vec\theta$ determines the MLE,
$\widehat{\vec\theta}$. The MLEs for the other parameters,
$\widehat{\vec\beta}$ and $\widehat{\sigma}$, are the corresponding
conditional estimates evaluated at $\widehat{\vec\theta}$.

### Determining $r^2_\theta$ and $\hat{\vec\beta}_\theta$ {#sec:betahat}

To determine $\tilde{\vec u}$ and $\widehat{\vec\beta}_\theta$
simultaneously we rearrange the terms in
([\[eq:PLSLMMP\]](#eq:PLSLMMP){reference-type="ref"
reference="eq:PLSLMMP"}) as $$\label{eq:PLSLMM1}
  \begin{bmatrix}
    \tilde{\vec u}\\
    \widehat{\vec\beta}_\theta
  \end{bmatrix}
  =\arg\min_{\vec u,\vec\beta}
  \left\|
    \begin{bmatrix}
      \vec y_{\text{obs}}\\
      \vec 0
    \end{bmatrix} -
    \begin{bmatrix}
      \vec Z\Lambda_\theta\vec P\trans & \vec X\\
      \vec P\trans &\vec 0
    \end{bmatrix}
    \begin{bmatrix}
      \vec P\vec u\\
      \vec\beta
    \end{bmatrix}
  \right\|^2 .$$ The PLS values, $\tilde{\vec u}$ and
$\widehat{\vec\beta}_\theta$, are the solutions to $$\label{eq:bigPLS}
  \begin{bmatrix}
    \vec P\left(\Lambda_\theta\trans\vec Z\trans\vec
      Z\Lambda_\theta+\vec I_q\right)\vec P\trans &
    \vec P\Lambda_\theta\trans\vec Z\trans\vec X\\
    \vec X\trans\vec Z\Lambda_\theta\vec P\trans & \vec X\trans\vec X
  \end{bmatrix}
  \begin{bmatrix}
    \vec P\tilde{\vec u}\\
    \widehat{\vec\beta}_\theta
  \end{bmatrix}=
  \begin{bmatrix}
    \vec P\Lambda_\theta\trans\vec Z\trans\vec y_{\text{obs}}\\
    \vec X\trans\vec y_{\text{obs}}
  \end{bmatrix}.$$ To evaluate these solutions we decompose the system
matrix as $$\label{eq:bigdecomp}
  \begin{bmatrix}
    \vec P\left(\Lambda_\theta\trans\vec Z\trans\vec
      Z\Lambda_\theta+\vec I_q\right)\vec P\trans &
    \vec P\Lambda_\theta\trans\vec Z\trans\vec X\\
    \vec X\trans\vec Z\Lambda_\theta\vec P\trans & \vec X\trans\vec X
  \end{bmatrix}
  =
  \begin{bmatrix}
    \vec L_\theta & \vec 0\\
    \vec R_{ZX}\trans & \vec R_X\trans
  \end{bmatrix}
  \begin{bmatrix}
    \vec L_\theta\trans & \vec R_{ZX}\\
    \vec 0 & \vec R_X
  \end{bmatrix}$$ where, as before, $\vec L_\theta$, the sparse Cholesky
factor, is the sparse lower triangular $q\times q$ matrix satisfying
([\[eq:sparseCholeskyP\]](#eq:sparseCholeskyP){reference-type="ref"
reference="eq:sparseCholeskyP"}). The other two matrices in
([\[eq:bigdecomp\]](#eq:bigdecomp){reference-type="ref"
reference="eq:bigdecomp"}): $\vec R_{ZX}$, which is a general
$q\times p$ matrix, and $\vec R_X$, which is an upper triangular
$p\times p$ matrix, satisfy $$\label{eq:RZXdef}
  \vec L_\theta\vec R_{ZX}=\vec P\Lambda_\theta\trans\vec Z\trans\vec X$$
and $$\label{eq:RXdef}
  \vec R_X\trans\vec R_X=\vec X\trans\vec X-\vec R_{ZX}\trans\vec R_{ZX}.$$

Those familiar with standard ways of writing a Cholesky decomposition as
either $\vec L\vec L\trans$ or $\vec R\trans\vec R$ ($\vec L$ is the
factor as it appears on the left and $\vec R$ is as it appears on the
right) will notice a notational inconsistency in
([\[eq:bigdecomp\]](#eq:bigdecomp){reference-type="ref"
reference="eq:bigdecomp"}). One Cholesky factor is defined as the lower
triangular fractor on the left and the other is defined as the upper
triangular factor on the right. It happens that in the Cholesky factor
of a dense positive-definite matrix is returned as the right factor,
whereas the sparse Cholesky factor is returned as the left factor.

One other technical point that should be addressed is whether $\vec
X\trans\vec X-\vec R_{ZX}\trans\vec R_{ZX}$ is positive definite. In
theory, if $\vec X$ has full column rank, so that $\vec X\trans\vec X$
is positive definite, then the downdated matrix, $\vec X\trans\vec
X-\vec R_{ZX}\trans\vec R_{ZX}$, must also be positive definite (see
Prob. [\[pr:th:downdate\]](#pr:th:downdate){reference-type="ref"
reference="pr:th:downdate"}). In practice, the downdated matrix can
become computationally singular in ill-conditioned problems, in which
case an error is reported.

The extended decomposition
([\[eq:bigdecomp\]](#eq:bigdecomp){reference-type="ref"
reference="eq:bigdecomp"}) not only provides for the evaluation of the
profiled deviance function, $\tilde{d}(\vec\theta)$,
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
  \bc Y\sim\mathcal{N}(\vec X\vec\beta,\sigma^2\vec I_n),$$ in which we
typically estimate $\sigma^2$ as $$\label{eq:21}
  \widehat{\sigma^2_R}=\frac{\|\vec y_{\text{obs}}-\vec X\widehat{\vec\beta}\|^2}{n-p}$$
even though the maximum likelihood estimate of $\sigma^2$ is
$$\label{eq:22}
  \widehat{\sigma^2_{L}}=\frac{\|\vec y_{\text{obs}}-\vec
    X\widehat{\vec\beta}\|^2}{n} .$$

The argument for preferring $\widehat{\sigma^2_R}$ to
$\widehat{\sigma^2_{L}}$ as an estimate of $\sigma^2$ is that the
numerator in both estimates is the sum of squared residuals at
$\widehat{\vec\beta}$ and, although the residual vector, $\vec
y_{\text{obs}}-\vec X\widehat{\vec\beta}$, is an $n$-dimensional vector,
it satisfies $p$ linearly independent constraints, $\vec
X\trans(\vec y_{\text{obs}}-\vec X\widehat{\vec\beta})=\vec 0$. That is,
the residual at $\widehat{\vec\beta}$ is the projection of the observed
response vector, $\vec y_{\text{obs}}$, into an $(n-p)$-dimensional
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
parameter estimates $\widehat{\vec\theta}_R$ and $\widehat{\sigma_R^2}$
for a linear mixed model have the property that they would specialize to
$\widehat{\sigma^2_R}$ from ([\[eq:21\]](#eq:21){reference-type="ref"
reference="eq:21"}) for a linear regression model, as seen in .

Although not usually derived in this way, the REML criterion (on the
deviance scale) can be expressed as $$\label{eq:23}
  d_R(\vec\theta,\sigma|\vec y_{\text{obs}})=-2\log
  \int_{\mathbb{R}^p}L(\vec\theta,\vec\beta,\sigma|\vec y_{\text{obs}})\,d\vec\beta .$$
The REML estimates $\widehat{\vec\theta}_R$ and $\widehat{\sigma_R^2}$
minimize $d_R(\vec\theta,\sigma|\vec y_{\text{obs}})$.

To evaluate this integral we form an expansion, similar to
([\[eq:PRSSwithL\]](#eq:PRSSwithL){reference-type="ref"
reference="eq:PRSSwithL"}), of $r^2_{\theta,\beta}$ about
$\widehat{\vec\beta}_\theta$ $$\label{eq:rsqbetathetaexp}
  r^2_{\theta,\beta}=r^2_\theta+\|\vec R_X(\vec\beta-\widehat{\vec\beta}_\theta)\|^2 .$$
In the same way that
([\[eq:PRSSwithL\]](#eq:PRSSwithL){reference-type="ref"
reference="eq:PRSSwithL"}) was used to simplify the integral in
([\[eq:hintegral\]](#eq:hintegral){reference-type="ref"
reference="eq:hintegral"}), we can derive $$\label{eq:betaintegral}
  \int_{\mathbb{R}^p}\frac{\exp\left(-\frac{r^2_{\theta,\beta}}{2\sigma^2}\right)}
  {(2\pi\sigma^2)^{n/2}|\vec L_\theta|} \,d\vec\beta=
  \frac{\exp\left(-\frac{r^2_\theta}{2\sigma^2}\right)}
  {(2\pi\sigma^2)^{(n-p)/2}|\vec L_\theta||\vec R_X|}$$ corresponding to
a REML criterion on the deviance scale of $$\label{eq:REMLdev}
  d_R(\vec\theta,\sigma|\vec y_{\text{obs}})=(n-p)\log(2\pi\sigma^2)+
  2\log\left(|\vec L_\theta||\vec R_X|\right)+\frac{r^2_\theta}{\sigma^2} .$$
Plugging in the conditional REML estimate,
$\widehat{\sigma^2}_R=r^2_\theta/(n-p)$, provides the profiled REML
criterion $$\label{eq:24}
  \tilde{d}_R(\vec\theta|\vec y_{\text{obs}})=
  2\log\left(|\vec L_\theta||\vec R_X|\right)+(n-p)
  \left[1+\log\left(\frac{2\pi r^2_\theta}{n-p}\right)\right].$$

The REML estimate of $\vec\theta$ is $$\label{eq:31}
  \widehat{\vec\theta}_R=\arg\min_{\vec\theta}\tilde{d}_R(\vec\theta|\vec y_{\text{obs}}) ,$$
and the REML estimate of $\sigma^2$ is the conditional REML estimate of
$\sigma^2$ at $\widehat{\vec\theta}_R$, $$\label{eq:REMLsigmasq}
  \widehat{\sigma^2_R}=r^2_{\widehat\theta_R}/(n-p) .$$ It is not
entirely clear how one would define a "REML estimate" of $\vec\beta$
because the REML criterion, $d_R(\vec\theta,\sigma|\vec
y)$, defined in ([\[eq:REMLdev\]](#eq:REMLdev){reference-type="ref"
reference="eq:REMLdev"}), does not depend on $\vec\beta$. However, it is
customary (and not unreasonable) to use
$\widehat{\vec\beta}_R=\widehat{\vec\beta}_{\widehat{\vec\theta}_R}$ as
the REML estimate of $\vec\beta$.

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
$\vec\theta=[1,0,1]\trans$

Without going into detail let us just point out that, because the
function is created within another function, , has access to the
structures in the model, , which internally is identified as . The
deviance evaluation uses the value of $\vec\theta$, called , and
information in three of the model slots: the random-effects structure,
called , the fixed-effects structure, called , and the response
structure, called . There is also a diagonal matrix of weights, called ,
available to this function but for our evaluation that matrix is the
identity and we will ignore it.

The function, with line numbers, is

Lines 3 to 9 create local versions of $\Lambda$ and $\vec L$ from the
information in the slot and the argument .

In lines 3 to 5 the argument is checked for correct mode and length and
whether it violates the lower bounds stored in the slot.

In lines 6 and 7 the current value of $\Lambda$ is created from a
template, stored in the slot, the value of and the index vector, which
maps elements of $\vec\theta$ to the non-zero elements of $\Lambda$.

As mentioned above, is the identity in this case so line 8 amounts to
assigning the matrix in the slot to the local variable . Then line 9
evaluates the sparse Cholesky factor, , by updating the template in the
slot. The optional argument, , to the method is the multiple of the
identity matrix to add to the of the second argument. This produces the
factor $\vec L_\theta$ defined in
([\[eq:sparseCholeskyP\]](#eq:sparseCholeskyP){reference-type="ref"
reference="eq:sparseCholeskyP"}).

Lines 10 to 14 produce $\vec Z\trans\vec y$, $\vec X\trans\vec y$ and
$\vec Z\trans\vec X$, taking into account possible weights, from the
(square root of the residual weights) slot, and/or an offset, from the
slot. To avoid confusion when incorporating weights, these values are
stored as , and .

Because the -based evaluation function, , must be able to handle
non-default values of arguments such as and to the function, we will
skip over some of the details and concentrate on the parts that
correspond to formulas in the previous section.

First we derive the matrix $\Lambda_\theta$ from $\vec\theta$. A
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
$\vec\theta$ (), the relative covariance factor, $\Lambda_\theta$ (),
the sparse Cholesky factor, $\vec
L_\theta$ (), the matrices $\vec R_{ZX}$ () and $\vec R_X$ (), the
conditional mode, $\tilde{\vec u}$ (), and the conditional estimate,
$\widehat{\vec\beta}_\theta$ (). The permutation represented by $\vec P$
is contained in the sparse Cholesky representation, .

Although the model matrices, $\vec X$ () and $\vec Z\trans$ (), and the
response vector, $\vec y_{\text{obs}}$ (), are available in the
environment, many of the products that involve only these fixed values
are precomputed and stored separately under the names
($\vec X\trans\vec X$), , and .

To provide easy access to the objects in the environment of we attach it
to the search path.

Please note that this is done here for illustration only. The practice
of attaching a list or a data frame or, less commonly, an environment in
an session is overused, somewhat dangerous (because of the potential of
forgetting to detach it later) and discouraged. The preferred practice
is to use the function to gain access by name to components of such
composite objects. For this section of code, however, using or would
quickly become very tedious and we use instead.

To update the matrix $\Lambda_\theta$ to a new value of $\vec\theta$ we
need to know which of the non-zeros in $\Lambda$ are updated from which
elements of $\vec\theta$. Recall that the dimension of $\vec\theta$ is
small (3, in this case) but $\Lambda$ is potentially large ($18\times18$
with $54$ non-zeros). The environment contains an integer vector that
maps the elements of to the non-zeros in .

Suppose we wish to recreate the evaluation of the profiled deviance at
the initial value of $\vec\theta=(1,0,1)$. We begin by updating
$\Lambda_\theta$ and forming the product
$\vec U\trans=\Lambda_\theta\trans\vec Z\trans$

The Cholesky factor object, , can be updated from without forming
$\vec U\trans\vec U+\vec I$ explicitly. The optional argument to the
method specifies a multiple of the identity to be added to
$\vec U\trans\vec U$

Then we evaluate and according to
([\[eq:RZXdef\]](#eq:RZXdef){reference-type="ref"
reference="eq:RZXdef"}) and
([\[eq:RXdef\]](#eq:RXdef){reference-type="ref" reference="eq:RXdef"})

Solving ([\[eq:bigPLS\]](#eq:bigPLS){reference-type="ref"
reference="eq:bigPLS"}) for $\tilde{\vec u}$ and
$\widehat{\vec\beta}_\theta$ is done in stages. Writing $\vec c_u$ and
$\vec c_\beta$ for the intermediate results that satisfy
$$\label{eq:stage1}
    \begin{bmatrix}
    \vec L_\theta & \vec 0\\
    \vec R_{ZX}\trans & \vec R_X\trans
  \end{bmatrix}
  \begin{bmatrix}
    \vec c_u\\\vec c_\beta
  \end{bmatrix}=
  \begin{bmatrix}
    \vec P\Lambda_\theta\trans\vec Z\trans\vec y_{\text{obs}}\\
    \vec X\trans\vec y_{\text{obs}} .
  \end{bmatrix}$$ we evaluate

The next set of equations to solve is $$\label{eq:stage2}
    \begin{bmatrix}
    \vec L_\theta\trans & \vec R_{ZX}\\
    \vec 0 & \vec R_X
  \end{bmatrix}
  \begin{bmatrix}
    \vec P\tilde{\vec u}\\
    \widehat{\vec\beta}_\theta
  \end{bmatrix}=
  \begin{bmatrix}
    \vec c_U\\\vec c_\beta
  \end{bmatrix}.$$

We can now create the conditional mean, , the penalized residual sum of
squares, , the logarithm of the square of the determinant of $\vec L$, ,
and the profiled deviance, which, fortuitously, equals the value shown
earlier.

to avoid later name clashes.

In terms of the calculations performed, these steps describe exactly the
evaluation of the profiled deviance in . The actual function for
evaluating the deviance, accessible as , ---FIXME--- is a slightly
modified version of what is shown above. However, the modifications are
only to avoid creating copies of potentially large objects and to allow
for cases where the model matrix, $\vec X$, is sparse. In practice,
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
distribution, $(\bc Y|\bc U=\vec u)$, is not required to be Gaussian but
does preserve some of the properties of the spherical Gaussian
conditional distribution $$(\bc Y|\bc U=\vec u)\sim\mathcal{N}
  (\vec Z\Lambda_\theta\vec u+\vec X\vec\beta,\sigma^2\vec I_n)$$ from
the linear mixed model. In particular, the components of $\bc Y$ are
*conditionally independent*, given $\bc U=\vec u$. Furthermore, $\vec u$
affects the distribution only through the conditional mean, which we
will continue to write as $\vec\mu$, and it affects the conditional mean
only through the linear predictor,
$\vec\gamma=\vec Z\Lambda_\theta\vec u+\vec X\vec\beta$.

Typically we do not have $\vec\mu=\vec\gamma$, however. The elements of
the linear predictor, $\vec\gamma$, can be positive or negative or zero.
Theoretically they can take on any value between $-\infty$ and $\infty$.
But many distributional forms used in GLMMs put constraints on the value
of the mean. For example, the mean of a Bernoulli random variable,
modeling a binary response, must be in the range $0<\mu<1$ and the mean
of a Poisson random variable, modeling a count, must be positive. To
achieve these constraints we write the conditional mean, $\vec\mu$, as a
transformation of the unbounded predictor, written $\vec\eta$. For
historical, and some theoretical, reasons the inverse of this
transformation is called the *link function*, written
$$\label{eq:linkfunction}
  \vec\eta=\vec g(\vec\mu) ,$$ and the transformation we want is called
the *inverse link*, written $\vec g^{-1}$.

Both $\vec g$ and $\vec g^{-1}$ are determined by scalar functions, $g$
and $g^{-1}$, respectively, applied to the individual components of the
vector argument. That is, $\vec\eta$ must be $n$-dimensional and the
vector-valued function $\vec\mu=\vec g^{-1}(\vec\eta)$ is defined by the
component functions $\mu_i=g^{-1}(\eta_i),\,i=1,\dots,n$. Among other
things, this means that the Jacobian matrix of the inverse link,
$\frac{d\vec\mu}{d\vec\eta}$, will be diagonal.

Because the link function, $\vec g$, and the inverse link, $\vec
g^{-1}$, are nonlinear functions (there would be no purpose in using a
linear link function) many people use the terms "generalized linear
mixed model" and "nonlinear mixed model" interchangeably. We reserve the
term "nonlinear mixed model" for the type of models used, for example,
in pharmacokinetics and pharmacodynamics, where the conditional
distribution is a spherical multivariate Gaussian
$$\label{eq:NLMMconddist}
  (\bc Y|\bc U=\vec u)\sim\mathcal{N}(\vec\mu, \sigma^2\vec I_n)$$ but
$\vec\mu$ depends nonlinearly on $\vec\gamma$. For NLMMs the length of
the linear predictor, $\vec\gamma$, is a multiple, $ns$, of $n$, the
length of $\vec\mu$.

Like the map from $\vec\eta$ to $\vec\mu$, the map from $\vec\gamma$ to
$\vec\mu$ has a "diagonal" property, which we now describe. If we use
$\vec\gamma$ to fill the columns of an $n\times
s$ matrix, $\Gamma$, then $\mu_i$ depends only on the $i$th row of
$\Gamma$. In fact, $\mu_i$ is determined by a nonlinear model function,
$f$, applied to the $i$ row of $\Gamma$. Writing
$\vec\mu=\vec f(\vec\gamma)$ based on the component function $f$, we see
that the Jacobian of $\vec f$, $\frac{d\vec\mu}{d\vec\gamma}$, will be
the vertical concatenation of $s$ diagonal $n\times n$ matrices.

Because we will allow for generalized nonlinear mixed models (GNLMMs),
in which the mapping from $\vec\gamma$ to $\vec\mu$ has the form
$$\label{eq:GammaEtaMu}
  \vec\gamma\;\rightarrow\;\vec\eta\;\rightarrow\;\vec\mu ,$$ we will
use ([\[eq:GammaEtaMu\]](#eq:GammaEtaMu){reference-type="ref"
reference="eq:GammaEtaMu"}) in our definitions.

### Determining the Conditional Mode, $\tilde{\vec u}$ {#sec:conditionalmode}

For all these types of mixed models, the conditional distribution,
$(\bc U|\bc Y=\vec y_{\text{obs}})$, is a continuous distribution for
which we can determine the unscaled conditional density, $h(\vec u)$. As
for linear mixed models, we define the conditional mode,
$\tilde{\vec u}$, as the value that maximizes the unscaled conditional
density.

Determining the conditional mode, $\tilde{\vec u}$, in a nonlinear mixed
model is a penalized nonlinear least squares (PNLS) problem
$$\label{eq:PNLSprob}
  \tilde{\vec u}=\arg\min_{\vec u}\|\vec y_{\text{obs}}-\vec\mu\|^2+\|\vec u\|^2$$
which we solve by adapting the iterative techniques, such as the
Gauss-Newton method [@bateswatts88:_nonlin Sect. 2.2.1], used for
nonlinear least squares. Starting at an initial value, $\vec
u^{(0)}$, (the bracketed superscript denotes the iteration number) with
conditional mean, $\vec\mu^{(0)}$, we determine an increment
$\vec\delta^{(1)}$ by solving the penalized linear least squares
problem, $$\label{eq:PNLSiter}
  \vec\delta^{(1)}=\arg\min_{\vec\delta}\left\|
    \begin{bmatrix}
      \vec y_{\text{obs}}-\vec\mu^{(0)}\\
      \vec 0-\vec u^{(0)}
    \end{bmatrix} -
    \begin{bmatrix}
      \vec U^{(0)}\\
      \vec I_q
    \end{bmatrix}\vec\delta\right\|^2$$ where $$\label{eq:Udef}
  \vec U^{(0)}=\left.\frac{d\vec\mu}{d\vec u}\right|_{\vec u^{(0)}} .$$
Naturally, we use the sparse Cholesky decomposition, $\vec
L_\theta^{(0)}$, satisfying $$\label{eq:sparseCholPNLS}
  \vec L_\theta^{(0)}\left(\vec L_\theta^{(0)}\right)=
  \vec P\left[\left(\vec U^{(0)}\right)\trans\vec U^{(0)}+\vec I_q\right]
  \vec P\trans$$ to determine this increment. The next iteration begins
at $$\label{eq:u1}
  \vec u^{(1)} = \vec u^{(0)} + k \vec\delta^{(1)}$$ where $k$ is the
step factor chosen, perhaps by step-halving [@bateswatts88:_nonlin
Sect. 2.2.1], to ensure that the penalized residual sum of squares
decreases at each iteration. Convergence is declared when the
orthogonality convergence criterion [@bateswatts88:_nonlin Sect. 2.2.3]
is below some pre-specified tolerance.

The *Laplace approximation* to the deviance is $$\label{eq:Laplace}
  d(\vec\theta,\vec\beta,\sigma|\vec y_{\text{obs}})\approx
  n\log(2\pi\sigma^2)+2\log|\vec L_{\theta,\beta}|+
  \frac{r^2_{\theta,\beta}}{\sigma^2},$$ where the Cholesky factor,
$\vec L_{\theta,\beta}$, and the penalized residual sum of squares,
$r^2_{\theta,\beta}$, are both evaluated at the conditional mode,
$\tilde{\vec u}$. The Cholesky factor depends on $\vec\theta$,
$\vec\beta$ and $\vec u$ for these models but typically the dependence
on $\vec\beta$ and $\vec u$ is weak.

Chapter Summary {#sec:lmmsummary}
---------------

The definitions and the computational results for maximum likelihood
estimation of the parameters in linear mixed models were summarized in .
A key computation is evaluation of the sparse Cholesky factor,
$\Lambda_\theta$, satisfying
([\[eq:sparseCholeskyP\]](#eq:sparseCholeskyP){reference-type="ref"
reference="eq:sparseCholeskyP"}), $$\vec L_\theta\vec L_\theta\trans
   = \vec P\left(\Lambda_\theta\trans\vec Z\trans\vec Z\Lambda_\theta
   +\vec I_q\right)\vec P\trans,$$ where $\vec P$ represents the
fill-reducing permutation determined during the symbolic phase of the
sparse Cholesky decomposition.

An extended decomposition
([\[eq:bigdecomp\]](#eq:bigdecomp){reference-type="ref"
reference="eq:bigdecomp"}) provides the $q\times p$ matrix $\vec R_{ZX}$
and the $p\times p$ upper triangular $\vec R_X$ that are used to
determine the conditional mode $\tilde{\vec u}_\theta$, the conditional
estimate $\widehat{\vec\beta}_\theta$, and the minimum penalized
residual sum of squares, $r^2_\theta$, from which the profiled deviance
$$\tilde{d}(\vec\theta|\vec y_{\text{obs}})
  =2\log|\vec L_\theta|+n\left[1 +
    \log\left(\frac{2 \pi r^2_\theta}{n}\right)\right]$$ or the profile
REML criterion $$\tilde{d}_R(\vec\theta|\vec y_{\text{obs}})=
  2\log\left(|\vec L_\theta||\vec R_X|\right)+(n-p)
  \left[1+\log\left(\frac{2\pi r^2_\theta}{n-p}\right)\right]$$ can be
evaluated and optimized (minimized) with respect to $\vec\theta$.

Exercises {#theory_exercises .unnumbered}
---------

Unlike the exercises in other chapters, these exercises establish
theoretical results, which do not always apply exactly to the
computational results.

[\[pr:th:posdef\]]{#pr:th:posdef label="pr:th:posdef"} Show that the
matrix $\vec A_\theta=\vec P\Lambda_\theta\trans\vec
  Z\trans \vec Z\Lambda_\theta\vec P\trans + \vec I_q$ is positive
definite. That is, $\vec b\trans\vec A\vec b>0,\,\forall\vec
  b\ne\vec 0$.

[\[pr:th:posdiag\]]{#pr:th:posdiag label="pr:th:posdiag"}

Show that $\Lambda_\theta$ can be defined to have non-negative diagonal
elements. (Hint: Show that the product $\Lambda_\theta\vec D$ where
$\vec D$ is a diagonal matrix with diagonal elements of $\pm1$ is also a
Cholesky factor. Thus the signs of the diagonal elements can be chosen
however we want.)

Use the result of
Prob. [\[pr:th:posdef\]](#pr:th:posdef){reference-type="ref"
reference="pr:th:posdef"} to show that the diagonal elements of
$\Lambda_\theta$ must be non-zero. (Hint: Suppose that the first zero on
the diagonal of $\Lambda_\theta$ is in the $i$th position. Show that
there is a solution $\vec x$ to $\Lambda_\theta\trans\vec x=\vec 0$ with
$x_i=1$ and $x_j=0,\,j=i+1,\dots,q$ and that this $\vec x$ contradicts
the positive definite condition.)

[\[pr:th:XtXposdef\]]{#pr:th:XtXposdef label="pr:th:XtXposdef"} Show
that if $\vec X$ has full column rank, which means that there does not
exist a $\vec\beta\ne\vec0$ for which $\vec X\vec\beta=\vec0$, then
$\vec X\trans\vec X$ is positive definite.

[\[pr:th:downdate\]]{#pr:th:downdate label="pr:th:downdate"} Show that
if $\vec X$ has full column rank then $$\begin{bmatrix}
      \vec Z\Lambda_\theta\vec P\trans & \vec X\\
      \vec P\trans &\vec 0
    \end{bmatrix}$$ also must have full column rank. (Hint: First show
that $\vec u$ must be zero in any vector $\begin{bmatrix}\vec
    u\\\vec\beta\end{bmatrix}$ satisfying $$\begin{bmatrix}
      \vec Z\Lambda_\theta\vec P\trans & \vec X\\
      \vec P\trans &\vec 0
    \end{bmatrix}
    \begin{bmatrix}\vec u\\\vec\beta\end{bmatrix}=\vec0 .$$ Use this
result and ([\[eq:bigdecomp\]](#eq:bigdecomp){reference-type="ref"
reference="eq:bigdecomp"}) to show that $$\begin{bmatrix}
      \vec L_\theta & \vec 0\\
      \vec R_{ZX}\trans & \vec R_X\trans
    \end{bmatrix}
    \begin{bmatrix}
      \vec L_\theta\trans & \vec R_{ZX}\\
      \vec 0 & \vec R_X
    \end{bmatrix}$$ is positive definite and, hence, $\vec R_X$ is
non-singular.)
