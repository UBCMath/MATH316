$$
\newcommand{\N}[1]{\left\|#1\right\|}
\newcommand{\abs}[1]{|#1|}
\newcommand{\mat}[1]{{\mathbf #1}}
\newcommand{\vect}[1]{\underline{#1}}
\newcommand{\njump}[1]{[|#1|]}
\newcommand{\pa}{\partial}
\newcommand{\De}{\Delta}
\newcommand{\ra}{\rightarrow}
\newcommand{\dst}{\displaystyle}
\newcommand{\la}{\lambda}
\newcommand{\al}{\alpha}
\newcommand{\rme}{{\rm e}}
\newcommand{\rmi}{{\rm i}}

\newcommand{\R}{\mathbb{R}}
\newcommand{\bke}[1]{\left ( #1 \right )}
\newcommand{\bkt}[1]{\left [ #1 \right ]}
\newcommand{\bket}[1]{\left \{ #1 \right \}}
\newcommand{\norm}[1]{\left \| #1 \right \|}
\newcommand{\bka}[1]{\left \langle #1 \right \rangle}
\renewcommand{\th}{\theta}
\newcommand{\om}{\omega}
\newcommand{\pd}{\partial}
\newcommand{\dis}{\displaystyle}
\newcommand{\ve}[1]{\mathbf{#1}}
\newcommand{\what}[1]{\widehat{#1}}
\newcommand{\tint}{\int\kern-.6em\int\kern-.6em\int}
$$

# Series Solutions To ODEs With Variable Coefficients

In this lecture we will introduce series methods for the solution of
variable coefficient ODE. We introduce the concepts of ordinary
points about which Taylor series solutions are obtained and singular
points about which more general solutions are required.

```{admonition} Key Concepts
Variable coefficient ODE, Series
Solutions, Ordinary Points and Taylor Series, Singular Points,
radius of convergence of power series.
```

## Series Solution Of Odes

### Power Series

$f(x)=a_0+a_1x+a_2x^2+\cdots +a_nx^n$ polynomial approximation.

Idea: Extend the polynomial to include infinitely many terms.

$$
\begin{eqnarray}
f(x) & = & a_0+a_1x+a_2x^2+\cdots +a_nx^n+\cdots\mbox{ Power Series}\nonumber\\
& = & \sum\limits_{n=0}^\infty a_nx^n
\end{eqnarray}
$$(ref0)

````{prf:example}
$\rme^x=1+\frac{x}{1!}+\frac{x^2}{2!}+\frac{x^3}{3!}+\cdots
+\frac{x^n}{n!}+\cdots =\sum\limits_{n=0}^\infty \frac{x^n}{n!}$.
````

More General Power Series:

$$
\begin{eqnarray}
f(x) & = & \sum\limits_{n=0}^\infty
a_n(x-x_0)^n=a_0+a_1(x-x_0)+a_2(x-x_0)^2+\cdots\
\end{eqnarray}
$$(ref1)

### Taylor Series

If we know all the derivatives of a function $f(x)$ at a single
point $x_0$ then we can solve for the coefficients of a power series
that represents the function at neighboring points $x$ as follows:

$$
\begin{eqnarray}
f(x) & = & \sum\limits_{n=0}^\infty a_n(x-x_0)^n=a_0+a_1(x-x_0)+a_2(x-x_0)^2+\cdots\nonumber\\
f^\prime (x) & = & a_1+2a_2(x-x_0)+3a_3(x-x_0)^2+\cdots +na_n(x-x_0)^n+\cdots\nonumber\\
& &\quad\Rightarrow f^\prime (x_0)=a_1\nonumber\\
f^{\prime\prime} (x) & = & 2a_2+3\cdot2a_3(x-x_0)+\cdots +n(n-1)a_n(x-x_0)^{n-2}+\cdots\nonumber\\
& &\quad\Rightarrow f^{\prime\prime} (x_0)=2a_2\nonumber\\
f^{(3)}(x) & = & 3!a_3+4\cdot 3\cdot 2(x-x_0)+\cdots +n(n-1)(n-2)a_n(x-x_0)^{n-3}+\cdots\nonumber\\
& &\quad\Rightarrow f^{(3)}(x_0)=3!a_3\nonumber\\
f^{(n)}(x_0) & = & n!a_n\Rightarrow a_n=\frac{f^{(n)}(x_0)}{n!}\nonumber\\
\mbox{Therefore } f(x) & = & \sum\limits_{n=0}^\infty
\frac{f^{(n)}(x_0)}{n!}(x-x_0)^n
\end{eqnarray}
$$(ref2)

Alternative Form of Taylor Series:

$$
\begin{equation}
f(x_0+h)=\sum\limits_{n=0}^\infty
\frac{f^{(n)}(x_0)}{n!}h^n
\end{equation}
$$(ref3)

````{prf:example} Taylor-Maclauren expansions of common functions.
$$
\begin{eqnarray}
\rme^x & = & \sum\limits_{n=0}^\infty \frac{x^n}{n!}\nonumber\\
\sin x  & = & \sum\limits_{n=0}^\infty (-1)^n\frac{x^{2n+1}}{(2n+1)!}\quad \sinh x=\sum\limits_{n=0}^\infty \frac{x^{2n+1}}{(2n+1)!}\\
\cos x & = & \sum\limits_{n=0}^\infty
(-1)^n\frac{x^{2n}}{(2n)!}\hspace{.5in} \cosh
x=\sum\limits_{n=0}^\infty \frac{x^{2n}}{(2n)!}\nonumber
\end{eqnarray}
$$(ref4)

$$
\begin{eqnarray}
\rme^{i\theta} & = & 1+i\theta +\frac{(i\theta )^2}{2!}+\frac{(i\theta )^3}{3!}+\cdots\nonumber\\
& = & \left( 1-\frac{\theta^2}{2!}+\frac{\theta^4}{4!}-\cdots\right) +i\left(\theta -\frac{\theta^3}{3!}+\cdots\right)\\
& = & \cos\theta +i \sin\theta\nonumber
\end{eqnarray}
$$(ref5)

$$
\begin{eqnarray}
\rme^{x} & = & 1+x +\frac{x^2}{2!}+\frac{x^3}{3!}+\cdots\nonumber\\
& = & \left( 1+\frac{x^2}{2!}+\frac{x^4}{4!}+\cdots\right) +\left(x +\frac{x^3}{3!}+\cdots\right)\nonumber\\
& = & \cosh x + \sinh x \nonumber
\end{eqnarray}
$$(ref6)

$$
\begin{eqnarray}
\rme^{-x} & = & \cosh x - \sinh x \nonumber\\
\cosh x & = & (\rme^{x}+\rme^{-x})/2\nonumber \\
\sinh x & = & (\rme^{x}-\rme^{-x})/2\nonumber
\end{eqnarray}
$$(ref7)
````

### Series Solution To A Constant Coefficient ODE

````{prf:example}
In this example we use power series to solve the linear ODE.

$$
\begin{equation}
y'+2y=0  \label{EG1}
\end{equation}
$$(EG1)

which we solved by integrating factor in the previous
lecture. Since the unknown solution $y(x)$ and all its derivatives
are defined implicitly by the ODE $y'=-2y$, let us look for a series
solution of the form: $y(x)=\sum\limits_{n=0}^\infty a_nx^n$.

$$
\begin{eqnarray}
y'&=&\sum\limits_{n=1}^\infty a_nnx^{n-1}\\
\mbox{Therefore }\quad y^{\prime}+2y&=&\sum\limits_{n=1}^\infty
a_nnx^{n-1}+\sum\limits_{n=0}^\infty 2a_nx^n=0\nonumber
\end{eqnarray}
$$(ref9)

In the first sum let

$$
\begin{eqnarray*}\begin{array}{ll}
m=n-1\quad &n=1\Rightarrow m=0\\
n=m+1 &\end{array}\end{eqnarray*}
$$(ref10)

$$
\begin{eqnarray}\begin{array}{rcl}
\mbox{Therefore }\sum\limits_{m=0}^\infty
a_{m+1}(m+1)x^m+\sum\limits_{n=0}^\infty 2a_nx^n &= &0\\
n\Leftrightarrow m:\quad\sum\limits_{m=0}^\infty\left\{
a_{m+1}(m+1)+2a_m\right\} x^m &= &0\end{array}
\end{eqnarray}
$$(ref11)

$$
\begin{eqnarray}\begin{array}{c}
a_{m+1}= - \frac{2}{(m+1)}a_m\\
\\
a_1=-2a_0,a_2= + \frac{2}{2}\frac{2}{1}a_0,a_3= - \frac{2}{3}\cdot\frac{2}{2}\cdot\frac{2}{1}a_0=(-1)^3\frac{2^3}{3!}a_0,\\
\\
\ldots ,a_m=(-1)^m\frac{2^m}{m!}a_0\\
\\
\mbox{Therefore }y(x)=a_0\sum\limits_{m=0}^\infty
\frac{(-2x)^m}{m!}=a_0\rme^{-2x}\end{array}
\end{eqnarray}
$$(ref12)
````

### Variable Coefficient ODEs: Ordinary And Singular Points

````{prf:example}
We consider the following Cauchy-Euler equation:

$$
\begin{equation}(x-1)y^{\prime\prime}+y^\prime
\label{eq:SingularPointExample} =0\end{equation}
$$(SingularPointExample)

whose solution is obtained as follows

$$
\begin{eqnarray}
y& = & (x-1)^r\Rightarrow r(r-1)+r=r^2=0\quad r=0,0.\nonumber\\
\mbox{So that } y(x) & = & A+B\ln |x-1| \mbox{ is the general solution} 
\label{eq:SingularPointExampleSol}
\end{eqnarray}
$$(SingularPointExampleSol)

___Method I:___

The first method we consider for obtaining a series solution to
{eq}`SingularPointExample` is to use the ODE to calculate all
the derivatives of $y(x)$ by direct differentiation and then to
substitute these derivatives into Taylor's formula. Although this
method can, in principle, be applied to any suitable ODE, we will
see that the computations can rapidly become tedious. However, this
method does highlight when the power series method will fail.

Assume that $y(0)$ and $y^{\prime }(0)$ are given, then

$$
\begin{eqnarray*}
y^{\prime \prime } &=&-\frac{y^{\prime }}{x-1}\Rightarrow y^{\prime
\prime
}(0)=y^{\prime }(0) \\
y^{(3)} &=&-\frac{y^{\prime \prime }}{x-1}+\frac{y^{\prime }}{\left(
x-1\right) ^{2}}\Longrightarrow y^{(3)}(0)=+y^{\prime \prime
}(0)+y^{\prime }(0)=2y^{\prime }(0)
\end{eqnarray*}
$$(ref15)

Substituting this into Taylor's formula $y(x)=y(0)+xy^{\prime }(0)+\frac{
x^{2}}{2}y^{\prime \prime }(0)+\frac{x^{3}}{3!}y^{(3)}(0)+\cdots $

$$
\begin{equation}
y(x)=y(0)+y^{\prime }(0)\left\{
x+\frac{x^{2}}{2}+\frac{x^{3}}{3}+\cdots \right\}
\label{eq:SingularPointExampleSol1}
\end{equation}
$$(SingularPointExampleSol1)

We observe that this process works for equation {eq}`SingularPointExample` 
using the expansion point $x_{0}=0,$ but will not work for $x_{0}=1,$ 
which is called a __singular point__. In fact, a power series expansion is possible for
all points $x_{0}\neq 1,$ which are called ordinary points.

___Method II:___

We now consider an alternative, and more convenient, method for
determining the coefficients of a power series solution to the ODE,
by substituting the power series into the ODE and determining
recursion formulae for these coefficients. Expand around the
ordinary point $x_0=0$:

$$
\begin{equation}
y(x)=\sum\limits_{n=0}^\infty c_nx^n,\quad y^\prime
=\sum\limits_{n=1}^\infty nc_nx^{n-1},\quad y^{\prime\prime}
=\sum\limits_{n=2}^\infty c_nn(n-1)x^{n-2}
\end{equation}
$$(ref17)

$$
\begin{eqnarray}
(x-1)\sum\limits_{n=2}^\infty c_nn(n-1)x^{n-2}+\sum\limits_{n=1}^\infty nc_nx^{n-1}=0\nonumber\\
-\sum\limits_{n=2}^\infty c_nn(n-1)x^{n-2}+\sum\limits_{n=2}^\infty c_n\left\{ n(n-1)+n\right\} x^{n-1}+c_1=0\\
m-1=n-2\Rightarrow m=n-1\quad n=2\Rightarrow m=1\quad n=m+1\nonumber\\
-c_2\cdot 2\cdot 1+c_1+\sum\limits_{m=2}^\infty \left[
-c_{m+1}(m+1)m+c_mm^2\right] x^{m-1}=0\nonumber
\end{eqnarray}
$$(ref18)

where $c_0$ and $c_1$ are arbitrary:

$$
\begin{eqnarray}
c_{m+1} & = & \frac{m}{m+1}c_m\quad m\geq 2\quad c_2=\frac{c_1}{2}\nonumber\\
c_3 & = & \frac{2}{3}c_2=\frac{c_1}{3}\quad
c_4=\frac{3}{4}c_3=\frac{c_1}{4}\ldots c_n=\frac{c_1}{n} \nonumber
\end{eqnarray}
$$(ref19)

Therefore

$$
\begin{equation}
y(x) = c_0+c_1\sum\limits_{n=1}^\infty \frac{x^n}{n}.
\label{eq:SingularPointExampleSol2}
\end{equation}
$$(SingularPointExampleSol2)
````

We observe that the solutions {eq}`SingularPointExampleSol1`
and {eq}`SingularPointExampleSol2` obtained by the two
different methods are identical.

````{admonition} Recall
:class: tip
$$
\begin{eqnarray}
\frac{1}{1-x} & = & 1+x+x^2+\cdots\quad\int\frac{1}{1-x}\, dx=-\ln |1-x|=x+\frac{x^2}{2}+\frac{x^3}{3}+\cdots\nonumber\\
y(x) & = & A+B\ln |x-1|
\end{eqnarray}
$$(ref21)
````

Thus the series solution is identical to the solution {eq}`SingularPointExampleSol` provided $|x|\le1$.
We note that the radius of convergence of convergence for the power series {eq}`SingularPointExampleSol1`
is $1$, which corresponds to the distance between the expansion point $x_{0}=0$
and the nearest singular point $x=1.$

### Power Series Solution of General Variable Coefficient Linear ODEs

Consider solving variable coefficient linear ODEs of the form:

$$
\begin{equation}
P(x)y^{\prime\prime}+Q(x)y^\prime +R(x)y=0\mbox{ Homogeneous Eq.}
\label{variablecoeff}
\end{equation}
$$(ref22)

Divide through by $P(x)$:

$$
\begin{equation}
Ly=y^{\prime\prime}+p(x)y^\prime +q(x)y=0\quad p(x)=Q/P, \quad q(x)
= R/P \label{LyoverPEQ}
\end{equation}
$$(LyoverPEQ)
In order to calculate the higher derivatives of $y(x)$ to substitute
into Taylor's formula, we rewrite {eq}`LyoverPEQ` as follows:
$y^{\prime \prime }=-p(x)y^{\prime }-q(x)y$.

 If $y(x_{0})$ and $y^{\prime }(x_{0})$ are given, then $y^{\prime
\prime }(x_{0})$ can be obtained directly from the ODE. Higher
derivatives of $y$ can, in turn, be obtained by differentiating the
ODE repeatedly. This process will be successful provided $p(x)$ and
$q(x)$ are infinitely differentiable at $x=x_{0}$. In this case
$p(x)$ and $q(x)$ are said to be analytic at $x_{0}$ and have Taylor
expansions of the form

$$
\begin{eqnarray*}\begin{array}{lrcl}
\mbox{i.e. } &p(x)=p_0+p_1(x-x_0)+\cdots &=
&\sum\limits_{k=0}^\infty p_k(x-x_0)^k\\
 &q(x)=q_0+q_1(x-x_0)+\cdots &= &\sum\limits_{k=0}^\infty
 q_k(x-x_0)^k\end{array}
\end{eqnarray*}
$$(ref24)

___Ordinary points:___

The expansion point $x_0$ is said to be an ordinary point of
{eq}`LyoverPEQ` if $p(x)=Q/P$ and $q(x)=R/P$ are analytic at
$x_0$. If $x_0$ is an ordinary point  it is possible to obtain power
series expansions of the solution $y(x)$ of the form:

$$
\begin{equation}
y(x)=\sum\limits_{n=0}^\infty c_n(x-x_0)^n. \label{PowerSeries1}
\end{equation}
$$(PowerSeries1)
 The idea is to substitute the expansion {eq}`PowerSeries1` into {eq}`LyoverPEQ` and solve for the unknown coefficients $c_n$
 in order to determine a solution.

````{prf:observation}
- If $P$, $Q$ and $R$ are polynomials then a point $x_0$ such that $P(x_0)\not= 0$ is an ordinary point.
- If $x_0=0$ is an ordinary point then we assume

$$
\begin{eqnarray}
y & = & \sum\limits_{n=0}^\infty c_nx^n, y_n^\prime =\sum\limits_{n=1}^\infty c_nnx^{n-1}, y_n^{\prime\prime}=\sum\limits_{n=2}^\infty c_nn(n-1)x^{n-2}\nonumber\\
0=Ly & = & \sum\limits_{n=2}^\infty c_nn(n-1)x^{n-2}+\left(\sum\limits_{n=0}^\infty p_n x^n\right) \sum\limits_{n=1}^\infty nc_nx^{n-1}\\
& &\quad +\, \left(\sum\limits_{n=0}^\infty q_n
x^n\right)\left(\sum\limits_{n=0}^\infty c_nx^n\right)\nonumber
\end{eqnarray}
$$(ref26)

$$
\begin{eqnarray}
& &\sum\limits_{m=0}^\infty\left\{ (m+2)(m+1)c_{m+2}+\big(
p_0(m+1)c_{m+1}+\cdots +p_mc_1\big)\right.\nonumber\\
& &\quad +\, \left.\left( q_0c_m+\cdots +q_mc_0\right)\right\} x^m=0
\end{eqnarray}
$$(ref27)
yields a non-degenerate recursion for the $c_m$. At an ordinary
point $x_0$ we can obtain two linearly independent solutions  of the
form {eq}`PowerSeries1`.
````

___Singular Points:___

If $p(x)$ or $q(x)$ are not analytic at $x_0$, then $x_0$ is said to
be a singular point of {eq}`LyoverPEQ`. For example if $P$, $Q$
and $R$ are polynomials and $P(x_0)=0$ and $Q(x_0)\not= 0$ or
$R(x_0)\not= 0$ then $x_0$ is a singular point. Or if
$p(x)=\sqrt(x)$ and $q(x)=2$, then $x_0=0$ is a singular point
because $p(x)$ is not differentiable at $x=0$.

The radius of convergence of {eq}`PowerSeries1` is at least as
large as the radius of convergence of each of the series expansions
for $p(x)=Q/P$ and $q(x)=R/P$, i.e., up to the closest singularity
to $x_0$.

````{prf:example} The Airy Equation
Consider the Airy equation, which arises in Quantum Mechanics:

$$
\begin{equation}
Ly = y^{\prime\prime} - xy=0 \label{eq:Airy}
\end{equation}
$$(ref28)

We observe that $x =0$ is an ordinary point.

$$
\begin{eqnarray}\begin{array}{c}
y=\sum\limits_{n=0}^\infty c_nx^n,\quad y^\prime =\sum\limits_{n=1}^\infty c_nnx^{n-1},\quad y^{\prime\prime}=\sum\limits_{n=2}^\infty c_nn(n-1)x^{n-2}\\
\sum\limits_{n=2}^\infty c_nn(n-1)x^{n-2}=\sum\limits_{n=0}^\infty c_nx^{n+1}\\
m+1=n-2\quad n=m+3\quad n=2\Rightarrow m=-1\\
c_22x^0+\sum\limits_{m=0}^\infty \big[ c_{m+3}(m+3)(m+2)-c_m\big] x^{m+1}=0\\
c_2=0\quad c_{m+3}=\frac{c_m}{(m+3)(m+2)}\quad
m=0,1,\ldots\end{array}
\end{eqnarray}
$$(ref29)


1. $c_0\ra c_3\ra c_6$.

$$
\begin{eqnarray}
c_3 & = & \frac{c_0}{3.2}, c_6=\frac{c_3}{6.5}=\frac{c_0}{6.5.3.2},c_9=\frac{c_0}{9.8.6.5.3.2}\nonumber\\
c_{3n} & = & \frac{c_0}{(3n)(3n-1)(3n-3)(3n-4)\ldots 9.8.6.5.3.2}\\
y_0(x) & = & 1+\frac{x^3}{3.2}+\frac{x^6}{6.5.3.2}+\cdots
+\frac{x^{3n}}{(3n)(3n-1)\ldots 3.2}+\ldots\nonumber
\end{eqnarray}
$$(ref30)

2. $c_1\ra c_4\ra c_7\ra c_10$.

$$
\begin{equation}
c_4=\frac{c_1}{4.3}\quad c_7=\frac{c_1}{7.6.4.3}\quad
c_{10}=\frac{c_1}{(10.9)(7.6)(4.3)}
\end{equation}
$$(ref31)

$$
\begin{eqnarray}
c_{3n+1} & = &\frac{c_1}{(3n+1)(3n)(3n-2)(3n-3)\ldots (7.6)(4.3)}\nonumber\\
y_1(x) & = & x+\frac{x^4}{4.3}+\frac{x^7}{7.6.4.3}+\cdots +\frac{x^{3n+1}}{(3n+1)(3n)\ldots 4.3}\\
y(x) & = & c_0y_0(x)+c_1y_1(x)\nonumber
\end{eqnarray}
$$(ref32)
````

___Radius of Convergence:___

$$
\begin{equation}
\lim_{m\ra\infty}\frac{c_{m+3}}{c_m}|x|^3
=\lim_{m\ra\infty}\frac{|x|^3}{(m+3)(m+2)}=0<1\quad\rho =\infty .
\end{equation}
$$(ref33)

See B\&D for expansion of Airy Solution about $x_0=1$:,
i.e. $y(x)=\sum a_n(x-1)^n$. It is useful to write $x=(x-1)+1$.

$$
\begin{equation}y^{\prime\prime}=(x-1)y+y\end{equation}
$$(ref34)

````{prf:example} The Hermite Equation
Consider the Hermite equation, which has application in Quantum
mechanics and numerical analysis:

$$
\begin{equation}
Ly=y^{\prime\prime}-2xy^\prime +\la y=0 \label{eq:Hermite}
\end{equation}
$$(ref35)

Since $x=0$ is an ordinary point let $\dst
y(x)=\sum\limits_{n=0}^\infty a_nx^n$ then

$$
\begin{equation}
Ly=\sum\limits_{n=2}^\infty
a_nn(n-1)x^{n-2}-2\sum\limits_{n=1}^\infty
a_nnx^n+\la\sum\limits_{n=0}^\infty a_nx^n=0.
\end{equation}
$$(ref36)

$$
\begin{eqnarray*}\begin{array}{lll}
m=n-2\ra n=m+2\quad &m\leftarrow n\quad &m\leftarrow n\\
n=2\Rightarrow m=0 & &\end{array}
\end{eqnarray*}
$$(ref37)

Therefore

$$
\begin{equation}
\sum\limits_{m=1}^\infty \big[ a_{m+2}(m+2)(m+1)-2a_mm+\la a_m\big]
x^m+\left[ a_2 2+\la a_0\right] x^0=0.
\end{equation}
$$(ref38)
$x^0:$

$$
\begin{equation} a_2=-\la a_0/2
\end{equation}
$$(ref39)
$x^m:$

$$
\begin{equation} a_{m+2}=\frac{(2m-\la )a_m}{(m+1)(m+2)} \quad m\geq 1
\end{equation}
$$(ref40)

$a_0$:

$$
\begin{eqnarray}
a_2 & = & -\frac{\la}{2}a_0, a_4=\frac{(4-\la )}{4.3}
a_2=\frac{(4-\la )(-\la )}{4.3.2}a_0,
a_6=\frac{(8-\la )(4-\la )(-\la )}{6.5.4.3.2}a_0\nonumber\\
&\phantom{=} & \quad a_{2k} =\frac{[4(k-1)-\la ][4(k-2)-\la ]\ldots
[-\la] }{(2k)!}a_0\\
y_0 & = & a_0\left[ 1-\frac{\la}{2}x^2+\frac{(\la -4)\la}{4!}
x^4+\frac{(8-\la )(4-\la )(-\la )}{6!}x^6+\cdots \right]\nonumber
\end{eqnarray}
$$(ref41)

$a_1$:

$$
\begin{eqnarray}
a_3 & = & \frac{(2-\la)}{3.2}a_1; a_5=\frac{(6-\la )}{5!}(2-\la )a_1; a_7=\frac{(10-\la )(6-\la )(2-\la )}{7!}a_1, \ldots\\
y_1 & = & a_1\left[ x+\frac{(2-\la )}{3!}x^3+\frac{(6-\la )(2-\la
)}{5!}x^5+\frac{(10-\la )(6-\la )(2-\la
)x^7}{7!}+\cdots\right]\nonumber
\end{eqnarray}
$$(ref42)

The general solution is of the form

$$
\begin{equation}
y(x)=Ay_0(x)+By_1(x)
\end{equation}
$$(ref43)
````

````{note}
(a) If $\la =2n$ then the recursion yields $a_{m+2}=0=a_{m+4}=\cdots $ for $m=n$.
Thus if $n$ is an even integer then the series solution $y_0$ will
terminate and become a polynomial of degree $n$.

In this case:

$$
\begin{eqnarray}
y_0(x) &=& a_0\left[ 1-nx^2+n(n-2)2^2\frac{x^4}{4!}-n(n-2)(n-4)\frac{2^3x^6}{6!}+\cdots\right.\nonumber\\
& &\quad \left. +\, (-1)^{n/2} n(n-2)\ldots 2.\big(2^{n/2}\big)
\frac{x^n}{n!}\right] .
\end{eqnarray}
$$(ref44)

On the other hand if $n$ is an odd integer then the series solution
$y_1(x)$ will terminate and become a polynomial of degree $n$. In
this case

$$
\begin{eqnarray}
y_1(x) & = & a_1\left[ x-2(n-1)\frac{x^3}{3!}+2^2(n-1)(n-3)\frac{x^5}{5!}\right.\nonumber\\
& &\quad -\, (n-1)(n-3)(n-5)2^3\frac{x^7}{7!}+\cdots\\
& &\quad \left. +\, (n-1)(n-3) \ldots
3.1(-2)^{\frac{(n-1)}{2}}\frac{x^n}{n!}\right]\nonumber
\end{eqnarray}
$$(ref45)

(b) For example in the special case $\la =4=2n$ then $n=2$.

$$
\begin{equation} y_0(x)=a_0[1-2x^2].\end{equation}
$$(ref46)
````
