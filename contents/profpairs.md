Examining likelihood contour projections {#chap:ProfPairs}
========================================

### Profile Pairs Plots {#sec:Profpairs}

A profiled deviance object, such as , not only provides information on
the sensitivity of the model fit to changes in parameters, it also tells
us how the parameters influence each other. When we re-fit the model
subject to a constraint such as, say, $\sigma_1=60$, we obtain the
conditional estimates for the other parameters --- $\sigma$ and
$\beta_0$ in this case. The conditional estimate of, say, $\sigma$ as a
function of $\sigma_1$ is called the *profile trace* of $\sigma$ on
$\sigma_1$. Plotting such traces provides valuable information on how
the parameters in the model are influenced by each other.

The *profile pairs* plot, obtained as

and shown in
Fig. [\[fig:fm01profpair\]](#fig:fm01profpair){reference-type="ref"
reference="fig:fm01profpair"} shows the profile traces along with
interpolated contours of the two-dimensional profiled deviance function.
The contours are chosen to correspond to the two-dimensional marginal
confidence regions at particular confidence levels.

Because this plot may be rather confusing at first we will explain what
is shown in each panel. To make it easier to refer to panels we assign
them $(x,y)$ coordinates, as in a Cartesian coordinate system. The
columns are numbered 1 to 3 from left to right and the rows are numbered
1 to 3 from bottom to top. Note that the rows are numbered from the
bottom to the top, like the y-axis of a graph, not from top to bottom,
like a matrix.

The diagonal panels show the ordering of the parameters: $\sigma_1$
first, then $\log(\sigma)$ then $\beta_0$. Panels above the diagonal are
in the original scale of the parameters. That is, the top-left panel,
which is the $(1,3)$ position, has $\sigma_1$ on the horizontal axis and
$\beta_0$ on the vertical axis.

In addition to the contour lines in this panel, there are two other
lines, which are the profile traces of $\sigma_1$ on $\beta_0$ and of
$\beta_0$ on $\sigma_1$. The profile trace of $\beta_0$ on $\sigma_1$ is
a straight horizontal line, indicating that the conditional estimate of
$\beta_0$, given a value of $\sigma_1$, is constant. Again, this is a
consequence of the simple model form and the balanced data set. The
other line in this panel, which is the profile trace of $\sigma_1$ on
$\beta_0$, is curved. That is, the conditional estimate of $\sigma_1$
given $\beta_0$ depends on $\beta_0$. As $\beta_0$ moves away from the
estimate, $\widehat{\beta}_0$, in either direction, the conditional
estimate of $\sigma_1$ increases.

We will refer to the two traces on a panel as the "horizontal trace" and
"vertical trace". They are not always perfectly horizontal and vertical
lines but the meaning should be clear from the panel because one trace
will always be more horizontal and the other will be more vertical. The
one that is more horizontal is the trace of the parameter on the y axis
as a function of the parameter on the horizontal axis, and vice versa.

The contours shown on the panel are interpolated from the profile zeta
function and the profile traces, in the manner described in
@bateswatts88:_nonlin [Chapter 6]. One characteristic of a profile
trace, which we can verify visually in this panel, is that the tangent
to a contour must be vertical where it intersects the horizontal trace
and horizontal where it intersects the vertical trace.

The $(2,3)$ panel shows $\beta_0$ versus $\log(\sigma)$. In this case
the traces actually are horizontal and vertical straight lines. That is,
the conditional estimate of $\beta_0$ doesn't depend on $\log(\sigma)$
and the conditional estimate of $\log(\sigma)$ doesn't depend on
$\beta_0$. Even in this case, however, the contour lines are not
concentric ellipses, because the deviance is not perfectly quadratic in
these parameters. That is, the zeta functions, $\zeta(\beta_0)$ and
$\zeta(\log(\sigma))$, are not linear.

The $(1,2)$ panel, showing $\log(\sigma)$ versus $\sigma_1$ shows
distortion along both axes and nonlinear patterns in both traces. When
$\sigma_1$ is close to zero the conditional estimate of $\log(\sigma)$
is larger than when $\sigma_1$ is large. In other words small values of
$\sigma_1$ inflate the estimate of $\log(\sigma)$ because the
variability that would be explained by the random effects gets
incorporated into the residual noise term.

Panels below the diagonal are on the $\zeta$ scale, which is why the
axes on each of these panels span the same range, approximately $-3$ to
$+3$, and the profile traces always cross at the origin. Thus the
$(3,1)$ panel shows $\zeta(\sigma_1)$ on the vertical axis versus
$\zeta(\beta_0)$ on the horizontal. These panels allow us to see
distortions from an elliptical shape due to nonlinearity of the traces,
separately from the one-dimensional distortions caused by a poor choice
of scale for the parameter. The $\zeta$ scales provide, in some sense,
the best possible set of single-parameter transformations for assessing
the contours. On the $\zeta$ scales the extent of a contour on the
horizontal axis is exactly the same as the extent on the vertical axis
and both are centered about zero.

Another way to think of this is that, if we would have profiled
$\sigma_1^2$ instead of $\sigma_1$, we would change all the panels in
the first column but the panels on the first row would remain the same.

Exercises {#exercises .unnumbered}
---------

Create a profile pairs plot for model fit in
Chap. [\[chap:ExamLMM\]](#chap:ExamLMM){reference-type="ref"
reference="chap:ExamLMM"} to the data. Does the shape of the deviance
contours in this model mirror those in
Fig. [\[fig:fm01profpair\]](#fig:fm01profpair){reference-type="ref"
reference="fig:fm01profpair"}?
