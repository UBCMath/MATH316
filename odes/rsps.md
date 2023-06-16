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

# Regular Singular Points

In this lecture we will define a Regular Singular Point about which
a Taylor series will not work. We will also introduce the concept of
the radius of convergence of the series and how it relates to the
coefficient of the highest derivative of the ODE.

```{admonition} Key Concepts
Ordinary Points and Regular Singular
Points, radius of convergence of power series.
```

## Radius Of Convergence and Nearest Singular Points

````{prf:example}
:label: example-odes-rsps-0
$$
(1+x^2)y^{\prime\prime}+2xy^\prime +4x^2y=0
$$

1. If we were given $y(0)=0$ and $y^\prime (0)=1$ then we would want a power series expansion of the form

$$
\begin{equation}
y=\sum\limits_{n=0}^\infty c_nx^n\quad\mbox{about }x_0=0.
\end{equation}
$$(ref-odes-rsps-0)

Roots of $1+x^2=0$ are $x=\pm i$, so we expect the radius of
convergence of the TS for $\frac{1}{1+x^2}$ to be $1$ since

$$
\begin{equation}
\frac{1}{1+x^2} = 1-x^2+x^4-\cdots\,
\lim\left|\frac{a_{n+2}}{a_n}\right|  = 1\quad \rho = 1
\end{equation}
$$(ref-odes-rsps-1)

2. If we were given $y(1)=1$, $y^\prime (1)=0$ then a power series expansion of
the form $\sum c_n(x-1)^n$ is required. In this case $\rho
=\sqrt{2}$.
````

````{prf:example}
:label: example-odes-rsps-1
$$
(x-1)(2x-1)y^{\prime\prime}+2xy^\prime -2y=0
$$

$x=0$ is an ordinary point. $x=1$ and $x=\frac{1}{2}$ are
singular points. One solution of this equation is

$$
\begin{equation}
y(x)=\frac{1}{x-1}=-(1+x+x^2+\cdots )\quad\rho =1.
\end{equation}
$$(ref-odes-rsps-2)

This Taylor Series solution about the ordinary point $x=0$ converges
beyond the singular point $x=\frac{1}{2}$.
````

````{prf:example}
:label: example-odes-rsps-2 
$$
(x^2-2x)y^{\prime\prime}+5(x-1)y^\prime +3y=0
$$
$$
y(1)=7 
y^\prime (1)=3
$$

$x=1$ is an ordinary point. $x=0$ is a singular point $\big[
(x-1)^2-1\big] y^{\prime\prime}+5(x-1)y^\prime +3y=0$.

Let $t=x-1$ so that $\frac{d}{dt}=\frac{d}{dx}$ and the equation is
transformed to

$$
\begin{eqnarray}\begin{array}{c}
(t^2-1){\dot{\dot y}}+5t \dot{y}+3y=0\nonumber\\
y=\sum\limits_{n=0}^\infty c_nt^n,\quad \dot{y} =\sum\limits_{n=1}^\infty c_nnt^{n-1},\quad {\dot{\dot y}} =\sum\limits_{n=2}^\infty c_nn(n-1)t^{n-2}\nonumber\\
\sum\limits_{n=2}^\infty n(n-1)c_nt^n-\sum\limits_{n=2}^\infty n(n-1)c_nt^{n-2}+5\sum\limits_{n=1}^\infty nc_nt^n+3\sum\limits_{n=0}^\infty c_nt^n=0 \\
m=n-2 \quad n=m+2\quad n=2=m=0\nonumber\\
\sum\limits_{m=2}^\infty \left[ -c_{m+2}(m+2)(m+1)+\left\{ m(m-1)+5m+3\right\} c_m\right] t^m\nonumber\\
-2c_2+3c_0+\left[ -c_3 3.2+5c_1+3c_1\right] t=0\end{array}
\end{eqnarray}
$$(ref-odes-rsps-3)

$$
\begin{eqnarray}\begin{array}{crcll}
t^0>   &c_2 &= &\frac{3}{2}c_0 &\\
\\
t^1 > &c_3 &= &\frac{8}{6}c_1=\frac{4}{3}c_1 &\\
\\
t^m> &c_{m+2} &= &\frac{c_m(m+1)(m+3)}{(m+1) (m+2)}\quad &m\geq
2.\end{array}
\end{eqnarray}
$$(ref-odes-rsps-4)

$c_0$:

$$
\begin{eqnarray}\begin{array}{lcll}
c_4 &= &\frac{5c_2}{4}=\frac{5}{4}\left(\frac{3}{2}\right) c_0,\quad &c_6=\frac{7}{6}c_4=\frac{7}{6}\frac{5}{4}\frac{3}{2}c_0\\
y_0(x) &= &\sum\limits_{n=0}^\infty \frac{357\ldots
(2n+1)}{246\ldots (2n)}(x-1)^{2n} &\end{array}
\end{eqnarray}
$$(ref-odes-rsps-5)

$c_1$:

$$
\begin{eqnarray}\begin{array}{lcll}
c_5 &= &\frac{6}{5}c_3=\frac{6}{5}\frac{4}{3}c_1\quad &c_{2n+1}=\frac{46\ldots 2n+2}{35\ldots 2n+1}c_1\nonumber\\
y_1(x) &= &\sum\limits_{n=0}^\infty \frac{46\ldots 2n+2}{35\ldots
2n+1}(x-1)^{2n+1} &\end{array}
\end{eqnarray}
$$(ref-odes-rsps-6)

$$
\begin{equation}
\lim\limits_{n\rightarrow\infty}\frac{c_{m+2}}{c_m}=\frac{m+3}{m+1}=1\quad
\rho =1
\end{equation}
$$(ref-odes-rsps-7)


$$
\begin{eqnarray}\begin{array}{c}
y(x)=c_0y_0(x)+c_1y_1(x)\\
y(1)=c_0=7\quad y^\prime (1)=c_1=3.
\end{array}\end{eqnarray}
$$(ref-odes-rsps-8)
````

### Singular Points

Consider

$$
\begin{equation}
P(x)y^{\prime\prime}+Q(x)y^\prime +R(x)y=0.
\end{equation}
$$(ref-odes-rsps-9)

If $P$, $Q$ and $R$ are polynomials without common factors then
singular points are points $x_0$ at which $P(x_0)=0$.

````{note}
At singular points the solution is not necessarily analytic.
````

````{prf:example}
:label: example-odes-rsps-3
$$
\begin{array}{l}x^2y^{\prime\prime}+xy^\prime =0\\
y=x^r\rightarrow r(r-1)+r=0\rightarrow y=c_1+c_2\ln x\end{array}
$$

The $x^2y^{\prime\prime}$ admits wild behaviour.
````

````{prf:example}
:label: example-odes-rsps-4
$$
\begin{array}{l}x^2y^{\prime\prime}-2y=0\\
y=x^r\rightarrow r(r-1)-2=0\quad r=2,-1\quad
y=c_1x^{2}+c_2x^{-1}\end{array}
$$

Again the $x^2y^{\prime\prime}$ admits wild behaviour.
````

````{prf:example}
:label: example-odes-rsps-5
$$
\begin{array}{l}x^2y^{\prime\prime}-2xy^\prime +2y=0\\
y=x^r\rightarrow r(r-1)-2r+2=0\quad r=1,2\quad y=c_1x+c_2x^2\end{array}
$$
In this case both solutions are analytic.
````

### Regular Singular Points - Polynomial Coefficients

Notice that all these cases are equidimensional equations for which
we can identify solutions of the form $x^r$ or $x^r\log x$. There is
a special class of singular points called regular singular points in
which the singularities are no worse than those in the
equidimensional equations.

$$
\begin{equation}
x^2 y^{\prime\prime}+\alpha {x}y^\prime +\beta y=0.
\end{equation}
$$(ref-odes-rsps-10)

If $P$, $Q$ and $R$ are polynomials and suppose $P(x_0)=0$ then
$x_0$ is a regular singular point if

$$
\begin{equation}
\lim\limits_{x\rightarrow x_0}(x-x_0)\frac{Q(x)}{P(x)}\quad\mbox{and}\quad
\lim\limits_{x\rightarrow x_0}(x-x_0)^2\frac{R(x)}{P(x)}\quad\mbox{are
finite.}
\end{equation}
$$(ref-odes-rsps-11)

$$
\begin{eqnarray}
\mbox{I.E.}\quad (x-x_0)\frac{Q(x)}{P(x)} & = &
p_0+p_1(x-x_0)+p_2(x-x_0)^2+\cdots\nonumber\\
& &\quad\quad \rightarrow \mbox{ singularity no worse than $\frac{1}{x-x_0}$}\\
(x-x_0)^2 \frac{R(x)}{P(x)} & = &
q_0+q_1 (x-x_0)+q_2(x-x_0)^2+\cdots\nonumber\\
& &\quad\quad \rightarrow \mbox{ singularity no worse than
$\frac{1}{(x-x_0)^2}$}\nonumber
\end{eqnarray}
$$(ref-odes-rsps-12)

````{prf:example}
:label: example-odes-rsps-6
$$
\begin{eqnarray}\begin{array}{c}
(1-x^2)y^{\prime\prime}-2xy^\prime +4y=0\\
P=1-x^2\quad P(\pm 1)=0\quad Q=-2x\quad R=4\\
\lim\limits_{x\rightarrow
1}(x-1)\frac{(-2x)}{(1-x)(1+x)}=1\quad\lim\limits_{x\rightarrow
1}(x-1)^2\frac{4}{(1+x)(1-x)}=0\end{array}
\end{eqnarray}
$$(ref-odes-rsps-13)

$x=1$ is a R.S.P. (similarly for $x=-1$).
````

````{prf:example}
:label: example-odes-rsps-7
$$
\begin{eqnarray}\begin{array}{c}
x^3y^{\prime\prime}-y=0\\
P(x)=x^3\quad Q=0\quad R=-1\\
\lim\limits_{x\rightarrow 0}x^2\left(\frac{-1}{x^3}\right)
=\infty\end{array}
\end{eqnarray}
$$(ref-odes-rsps-14)

Thus $x=0$ is an irregular singular point. Actually $y \sim
x^{3/4}\{\rm\ e\}^{\pm 2/x^{1/2}}$ as ${x\rightarrow 0+}$ which is much wilder
than the simple power law $x^r$ or $x^r\log x$.
````

````{note}
Any singular point that is not a regular singular
point is called an irregular singular point.
````

````{prf:example}
:label: example-odes-rsps-8
$$
2(x-2)^2xy^{\prime\prime}+3xy^\prime +(x-2)y=0
$$ 

Singular points at $x=0,2$. $x=0$ is a regular singular point. $x=2$ is an irregular singular point.
````

### More General Definition Of A Regular Singular Point:

If $P$, $Q$, and $R$ are not limited to polynomials then consider

$$
\begin{eqnarray}\begin{array}{lrcl}
&P(x)y^{\prime\prime}+Q(x)y^\prime +R(x)y &= &0\\
\mbox{or}&&&\\
&x^2y^{\prime\prime}+x\left(\frac{xQ}{P}\right)y^\prime
+\left(\frac{x^2R}{P}\right) y&= &0\end{array}
\end{eqnarray}
$$(ref-odes-rsps-15)

$x=0$ is a regular singular point if
$\displaystyle\left(\frac{xQ}{P}\right)$ and
$\displaystyle\left(\frac{x^2R}{P}\right)$ are analytic at $x=0$, i.e.,

$$
\begin{equation}
\frac{xQ}{P}=p(x)=p_0+p_1x+\cdots\quad\mbox{and}\quad
\frac{x^2R}{P}=q(x)=q_0+q_1x+\cdots .
\end{equation}
$$(ref-odes-rsps-16)

In this case

$$
\begin{equation}
Ly=x^2y^{\prime\prime}+xp_0y^\prime +q_0y+ \overbrace{x\left\{
p_1xy^\prime +q_1y+\cdots\right\}}^{\mbox{small as }x\rightarrow 0} = 0.
\end{equation}
$$(ref-odes-rsps-17)

Then as $x\rightarrow 0$ $x^2y^{\prime\prime}+xp_0y^\prime
+q_0y\approx 0$ which is an Euler Equation which has solutions of
the form $y=x^r$.

Thus about a regular singular point we look for solutions of the
form $\displaystyle y=x^r\sum\limits_{n=0}^\infty
a_nx^n=\sum\limits_{n=0}^\infty a_nx^{n+r}$.

Our task is to determine:

1. $r$,
2. the coefficients $a_n$,
3. the radius of convergence.

````{prf:example}
:label: example-odes-rsps-9
$$
x^2y^{\prime\prime}+2(\{\rm\ e\}^x-1)y^\prime
+\{\rm\ e\}^{-x}\cos xy=0,\quad P=x^2$, $Q=2(\{\rm\ e\}^x-1)$, $R=\{\rm\ e\}^{-x}\cos
x
$$

$x=0$ is a singular point.

$$
\begin{eqnarray}
\lim\limits_{x\rightarrow 0}\frac{xQ}{P} & = & \lim\limits_{k\rightarrow 0}x
\frac{2(\{\rm\ e\}^x-1)}{x^2} = \lim\limits_{x\rightarrow 0}\frac{2(\{\rm\ e\}^x-1)}{x}
\stackrel{\frac{0}{0}}{=}\lim\limits_{x\rightarrow 0}\frac{2\{\rm\ e\}^x}{1}=2\mbox{ L'Hopital}\nonumber\\
\lim\limits_{x\rightarrow 0}\frac{x^2R}{P} & = & \lim\limits_{x\rightarrow
0}x^2\frac{\{\rm\ e\}^{-x}\cos x}{x^2}=1<\infty .
\end{eqnarray}
$$(ref-odes-rsps-18)

Since the quotient functions $p=xQ/P$ and $q=x^2R/P$ have Taylor
Expansions about $x=0$, $x=0$ is a regular singular point.
````
