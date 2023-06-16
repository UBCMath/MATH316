$$
\newcommand{\N}[1]{\left\|#1\right\|}
\newcommand{\abs}[1]{|#1|}
\newcommand{\mat}[1]{{\mathbf #1}}
\newcommand{\vect}[1]{\underline{#1}}
\newcommand{\njump}[1]{[|#1|]}
\newcommand{\bke}[1]{\left ( #1 \right )}
\newcommand{\bkt}[1]{\left [ #1 \right ]}
\newcommand{\bket}[1]{\left \{ #1 \right \}}
\newcommand{\norm}[1]{\left \| #1 \right \|}
\newcommand{\bka}[1]{\left \langle #1 \right \rangle}
\newcommand{\ve}[1]{\mathbf{#1}}
\newcommand{\what}[1]{\widehat{#1}}
$$

# Frobenius Series About Regular Singular Points

In this lecture we will summarize the classification of expansion
points $x_0$ for series as ordinary points for which Taylor Series
approximations are appropriate, regular singular points for which
Frobenius series expansions will work, and irregular singular points
for which neither power series expansions work. We also discuss the
radius of convergence of series expansions of ODE, which is at least
as large as the minimum distance from $x_0$ to the nearest other
singularity in the complex plane.

```{admonition} Key Concepts
Series Solutions; Ordinary Points and
Taylor Series; Regular Singular Points and Frobenius Series;
Irregular Singular Points; radii of convergence of power series
solutions of ODE.
```

## Fobenius Series Expansions About Regular Singular Points

### Series Expansion Summary:

Consider

$$
\begin{equation}
P(x)y^{\prime\prime}+Q(x)y^\prime +R(x)y=0
\end{equation}
$$(ref-odes-frobenius-0)
Divide by $P(x)$:

$$
\begin{equation}
L y = y^{\prime\prime}+
   \frac{Q(x)}{P(x)}y^\prime +
   \frac{R(x)}{P(x)}y=0, \label{eq:OPLequation}
\end{equation}
$$(OPLequation)

and define the functions $p(x)=\frac{Q(x)}{P(x)}$ and
$q(x)=\frac{R(x)}{P(x)}$. We observe that {eq}`OPLequation` can
be used to generate $y''(x_0)$ and all higher derivatives by
differentiating the ODE repeatedly. This can be used to generate a
Taylor series expansion, provided $P(x_0)\ne 0$, or more genreally
if $p(x)$ and $q(x)$ are analytic at $x_0$, which leads to the
following definition.

___Ordinary Points:___

$x_0$ is an ordinary point if $p(x)$ and $q(x)$ are analytic at
$x_0$ i.e.,

$$
\begin{eqnarray}
p(x) & = & p_0+p_1(x-x_0)+\cdots\nonumber\\
q(x) & = & q_0+q_1(x-x_0)+\cdots \label{eq:OPexpansionofpandq}.
\end{eqnarray}
$$(OPexpansionofpandq)

About an ordinary point $x_0$ we can obtain 2 linearly independent
solutions of the form

$$
\begin{equation}
y(x)=\sum\limits_{n=0}^\infty a_n(x-x_0)^n
\end{equation}
$$(ref-odes-frobenius-3)

whose radii of convergence are at least as large as those of $p$ and
$q$ in {eq}`OPexpansionofpandq` - i.e., the circle of
convergence extends at least as far as the singularity closest to
$x_0$ in the complex plane.

___Singular Points:___

If $P(x_0)=0$ then $p(x)$ and $q(x)$ may
fail to be analytic in which case $x_0$ is a singular point.

Regular Singular Points - an example:

To motivate how to proceed near singular points let us
consider the following example.

````{prf:example}
:label: rsp-example

$$
\begin{eqnarray}
Ly & = & 2x^2y^{\prime\prime}-xy^\prime +(1+x)y=0\quad x=0\mbox{ is a RSP.}\nonumber\\
\end{eqnarray}
$$(ref-odes-frobenius-4)

We observe that this equation is ``almost Equidimensional" but for the additional
$xy$ term. Let us rearrange this equation as follows

$$
\begin{equation}
L_0y = 2x^2y^{\prime\prime}-xy^\prime +y=-xy \label{recursion}
\end{equation}
$$(recursion)

As $x\rightarrow 0$ we observe that $xy\ll y$, so that the
term on the right side of {eq}`recursion` is much less than those
on the left. Thus (throwing away the small term for the moment) we
solve the homogeneous equidimensional equation $L_0y=0$

$$
\begin{equation}
L_0 y = 2x^2y^{\prime\prime}-xy^\prime +y=0,\  \mbox{whose general
solution is } y\left( x\right) =C_{1}x+C_{2}\sqrt{x}
\end{equation}
$$(ref-odes-frobenius-6)

Now take this solution and substitute it on the right side of
{eq}`recursion`

$$
\begin{equation}
L_0 y=2x^{2}y^{\prime \prime }-xy^{\prime }+y=-x\left(
C_{1}x+C_{2}\sqrt{x}\right),\  \mbox{whose general solution is }
y\left( x\right) =C_{1}\left( x-\frac{1}{3} x^{2}\right)
+C_{2}\sqrt{x}\left( 1-x\right)
\end{equation}
$$(ref-odes-frobenius-7)

Substituting this term again on the right side of {eq}`recursion` we obtain

$$
\begin{eqnarray}
L_0 y & = & 2x^2y^{\prime\prime}-xy^\prime +y=-\left(C_{1}\left(
x-\frac{1}{3} x^{2}\right)+C_{2}\sqrt{x}\left( 1-x\right)\right)\nonumber\\
\mbox{whose general solution is } y\left( x\right) &=&C_{1}\left(
x-\frac{1}{3}x^{2}+\frac{1}{30}x^{3}\right) +C_{2}\sqrt{x}\left(
1-x+\frac{1}{6}x^{2}\right)
\end{eqnarray}
$$(ref-odes-frobenius-8)

We observe that this procedure is generates a series
expansion of the form:

$$
\begin{equation}
y(x) = x^r \sum\limits_{n=0}^\infty a_n x^{n} \label{Frobenius}
\end{equation}
$$(Frobenius)

where $r$ is one of the roots of the indicial equation
associated with the equidimensional part of the ODE.
````

````{prf:observation}
1. Singular points that can be treated by this procedure are
known as Regular Singular Points

2. The form of the series {eq}`Frobenius` that is suitable for
determining the behavior of the solutions to ODE around such regular
singular points are called Frobenius Series.

3. Rather than proceed with this recursive approach, which can
rapidly become complicated, we will adopt a procedure in which we
substitute the series of the form {eq}`Frobenius` directly into
the ODE and solve for the unknown coefficients. This procedure is
illustrated in {prf:ref}`frobenius-example` below, which yields the same solution.
````

___Regular Singular points for general second order linear ODE___

Motivated by the previous example we consider the more general second
order linear ODE. In this case we consider $(x-x_0)^2 L y$

$$
\begin{equation}
(x-x_0)^2 L y= (x-x_0)^2y^{\prime\prime} +(x-x_0)\left\{ (x-x_0)
   \frac{Q(x)}{P(x)}\right\} y^\prime
   +\left\{(x-x_0)^2\frac{R(x)}{P(x)}\right\}y=0 \label{eq:RSPLequation}
\end{equation}
$$(RSPLequation)

and define coefficient functions $p(x)=\left\{ (x-x_0) \frac{Q(x)}{P(x)}\right\}$
and $q(x)=\left\{(x-x_0)^2\frac{R(x)}{P(x)}\right\}$. A point $x_0$ is a
regular singular point if $p(x)$ and $q(x)$ are both analytic at $x_0$, i.e.,

$$
\begin{eqnarray}
p(x)&=&(x-x_0)\frac{Q(x)}{P(x)}= p_0+p_1(x-x_0)+\cdots \nonumber\\
q(x)&=&(x-x_0)^2\frac{R(x)}{P(x)}= q_0+q_1(x-x_0)+\cdots
\label{eq:RSPexpansionofpandq}
\end{eqnarray}
$$(RSPexpansionofpandq)

Now substituting the expansions defined in
{eq}`RSPexpansionofpandq` into {eq}`RSPLequation` and
collecting powers of $(x-x_0)$, we obtain:

$$
\begin{equation}
(x-x_0)^2 Ly=(x-x_0)^2y^{\prime\prime}+(x-x_0)p_0y^\prime +q_0y+
\overbrace{(x-x_0)\left\{ p_1 (x-x_0) y^\prime +q_1y+p_2 (x-x_0)^2
y'+q_2(x-x_0)y+\cdots\right\}}^{\mbox{small as }x\rightarrow x_0} = 0
\label{eq:expandedLoperatorforRSP}
\end{equation}
$$(expandedLoperatorforRSP)

Thus, in the limit as $x \rightarrow x_0$ the operator $(x-x_0)^2 Ly
\approx L_0 y$, where

$$
\begin{equation}
L_0 y=(x-x_0)^2y^{\prime\prime}+(x-x_0)p_0y^\prime +q_0y=0
\label{eq:L0equation}
\end{equation}
$$(L0equation)

This implies that close to the expansion point $x_0$, the operator $L y$ has singularities no worse than
the Euler Equation {eq}`L0equation`. Since the higher order
terms (designated as small in {eq}`expandedLoperatorforRSP`)
cannot introduce more singular terms but rather corrections that are
higher powers in $(x-x_0)$, we are motivated to look for solutions
of the form

$$
\begin{equation}
y(x)=(x-x_0)^r\sum\limits_{n=0}^\infty
c_n(x-x_0)^n=\sum\limits_{n=0}^\infty c_n(x-x_0)^{n+r},
\end{equation}
$$(ref-odes-frobenius-14)

which is known as a Frobenius Series.

___Radius of convergence:___

The radius of convergence of this series is greater than or equal
to the distance between $x_0$ and the nearest singular point in the
complex plane.

___Irregular Singular Points:___

If $x_0$ is a singular point and
$p(x)=(x-x_0)\frac{Q(x)}{P(x)}$ and $(x-x_0)^2\frac{R(x)}{P(x)}$ are
not analytic, then $x_0$ is called an irregular singular point.

````{prf:example}
:label: example-odes-frobenius-0
$x=0$ is an irregular point of the second order equation 

$$
Ly=x^2 y''+(1+3x)y'+y=0
$$  

In this case the leading behavior of $y(x)$ as $x \rightarrow 0$ is $y(x)\approx c x e^{1/x}$,
which could not be captured by a Frobenius expansion.
````

````{prf:example}
:label: example-odes-frobenius-1
$x=0$ is an irregular point of the first order equation 

$$
Ly=x^2 y'+y=0
$$

The solution of this first order linear equation can be obtained by means of an integrating factor 
$F=e^{-1/x}$, which yields the solution $y(x)= c e^{1/x}$, which could not be captured by a 
Frobenius expansion about $x_0=0$.
````

### Frobenius Series Expansion

````{prf:example}
:label: frobenius-example

We revisit {prf:ref}`rsp-example` by using a Frobenius series to solve the
equations directly.

$$
\begin{eqnarray}
Ly & = & 2x^2y^{\prime\prime}-xy^\prime +(1+x)y=0\quad x=0\mbox{ is a RSP.}\nonumber\\
y & = & \sum\limits_{n=0}^\infty a_nx^{n+r}
\end{eqnarray}
$$(ref-odes-frobenius-15)

$$
\begin{eqnarray}
Ly & = & 2x^2\sum\limits_{n=0}^\infty a_n(n+r)(n+r-1)x^{n+r-2}
  -x\sum\limits_{n=0}^\infty a_n(n+r)x^{n+r-1}\nonumber\\
& &\quad\quad +\, (1+x)\sum\limits_{n=0}^\infty a_nx^{n+r} = 0\nonumber\\
& &\sum\limits_{n=0}^\infty a_n\left\{ 2(n+r)(n+r-1)-(n+r)+1\right\}
   x^{n+r}\nonumber\\
& &\quad\quad +\, \sum\limits_{n=0}^\infty a_nx^{n+r+1}=0\\
& & \quad\quad m=n+1\quad n=0\rightarrow m=1\nonumber\\
& & \quad\quad n=m-1\nonumber\\
\mbox{Therefore }& & a_0\left\{ 2r(r-1)-r+1\right\} x^r +
   \sum\limits_{n=1}^\infty \left[ a_n\left\{ 2(n+r)(n+r-1)\right.\right.\nonumber\\
& &\quad\quad -\, \left.\left. (n+r)+1\right\} +a_{n-1}\right]
x^{n+r} =0.\nonumber
\end{eqnarray}
$$(ref-odes-frobenius-16)

$x^r>$ Indicial Equation: $2r^2-3r+1=(2r-1)(r-1)=0\quad
r=\frac{1}{2}, \quad r=1$.

$a_0$ arbitrary

Recursion

$$
\begin{equation} a_n=\frac{-a_{n-1}}{(2n+2r-3)(n+r)+1}\end{equation}
$$(ref-odes-frobenius-17)

Let $r=1/2$:

$$
\begin{eqnarray}
a_n & = & \frac{-a_{n-1}}{(2n-2)(n+1/2)+1}
   =\frac{-a_{n-1}}{(n-1)(2n+1)+1}=\frac{-a_{n-1}}{n(2n-1)}\nonumber\\
n=1:\quad a_1 & = & \frac{-a_0}{1};\quad n=2: a_2=\frac{-a_1}{2.3}=\frac{+a_0}{2.3}\nonumber\\
a_3 & = & \frac{-a_2}{3.5}=\frac{-a_0}{1.(2.3)(3.5)};\quad
   a_4=\frac{-a_3}{4.7}=\frac{+a_0}{1(2.3)(3.5)(4.7)}\\
a_n & = & \frac{(-1)^na_0}{n!1.3.5.(2n-1)}=\frac{(-1)^n2^{(n-1)}a_0}{n(2n-1)!}\nonumber\\
y_1(x) &= & x^{1/2}\sum\limits_{n=0}^\infty
\frac{(-1)^n2^{(n-1)}}{n(2n-1)!}x^n\nonumber
\end{eqnarray}
$$(ref-odes-frobenius-18)

$r=1$:

$$
\begin{eqnarray}
a_n & = & \frac{-a_{n-1}}{(2n-1)(n+1)+1}=\frac{-a_{n-1}}{(2n+1)n}\nonumber\\
a_1 & = & \frac{-a_0}{3.1},\quad
a_2=\frac{-a_1}{5.2}=\frac{+a_0}{(1.3)(2.5)};\quad
   a_3=\frac{-a_2}{3.7}=\frac{-a_0}{(1.3)(2.5)(3.7)}\nonumber\\
a_n & = & \frac{(-1)^na_0}{n!3.5.7. (2n+1)}=
   \frac{(-1)^n2^na_0}{(2n+1)!}\\
y_2(x) & = & x\sum\limits_{n=0}^\infty
\frac{(-1)^n2^n}{(2n+1)!}x^n\nonumber
\end{eqnarray}
$$(ref-odes-frobenius-19)

__General Solution:__ $y(x)=c_1y_1(x)+c_2y_2(x)$
Radius of Convergence $\infty$.
````
