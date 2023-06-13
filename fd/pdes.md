$$
\newenvironment{remark}{\begin{Remark}\rm}{\end{Remark}}
\newenvironment{exercise}{\begin{Exercise}\rm}{\end{Exercise}}
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

# Solving the Heat, Laplace and Wave Equations Using Finite Difference Methods

In this lecture we introduce the finite difference method that is
widely used for approximating PDEs using the computer. We use the
definition of the derivative and Taylor series to derive finite
difference approximations to the first and second derivatives of a
function. We then use these finite difference quotients to
approximate the derivatives in the heat equation and to derive a
finite difference approximation to the heat equation. Similarly, the
technique is applied to the wave equation and Laplace's Equation.
The technique is illustrated using EXCEL spreadsheets.

```{admonition} Key Concepts
Finite Difference Approximations to
derivatives, The Finite Difference Method, The Heat Equation, The
Wave Equation, Laplace's Equation.
```

## Finite Difference Methods

(fd-intro)=

### Approximating the Derivatives of a Function by Finite Differences

Recall that the derivative of a function was defined by taking the
limit of a difference quotient:

$$
\begin{eqnarray}
f'(x)=\lim\limits_{\De x\ra 0}\frac{f(x+\De x)-f(x)}{\De x}
\label{Derivative}
\end{eqnarray}
$$(Derivative)

Now to use the computer to solve differential equations we go in the
opposite direction - we replace derivatives by appropriate
difference quotients. If we assume that the function can be
differentiated many times then Taylor's Theorem is a very useful
device to determine the appropriate difference quotient to use. For
example, consider

$$
\begin{eqnarray}
f(x+\De x)=f(x)+\De xf'(x)+\frac{\De x^2}{2!}f''(x)+\frac{\De
x^3}{3!}f^{(3)}(x)+\frac{\De x^4}{4!}f^{(4)}(x)+\ldots
\label{Taylor}
\end{eqnarray}
$$(Taylor)

Re-arranging terms in {eq}`Taylor` and dividing by $\De x$ we
obtain

$$
\begin{eqnarray*}
\frac{f(x+\De x)-f(x)}{\De x}=f'(x)+\frac{\De x}{2}f''(x)+\frac{\De
x^2}{3!}f^{(3)}(x)+\ldots .
\end{eqnarray*}
$$(ref2)

If we take the limit $\De x\ra 0$ then we recover
{eq}`Derivative`. But for our purposes it is more useful to retain
the approximation

$$
\begin{eqnarray}
\frac{f(x+\De x)-f(x)}{\De x} & = & f'(x)+\frac{\De x}{2}f''(\xi )\\
& = & f'(x)+O(\De x)\label{ForwardDiff} .
\end{eqnarray}
$$(ForwardDiff)

We retain the term $\dis\frac{\De x}{2}f''(\xi )$ in
{eq}`ForwardDiff` as a measure of the error involved when we
approximate $f'(x)$ by the difference quotient $\big( f(x+\De
x)-f(x)\big) /\De x$. Notice that this error depends on how large
$f''$ is in the interval $[x,x+\De x]$ (i.e. on the smoothness of
$f$) and on the size of $\De x$. Since we like to focus on that part
of the error we can control we say that the error term is of the
order $\De x$ -- denoted by $O(\De x)$. Technically a term or
function $E(\De x)$ is $O(\De x)$ if

$$
\begin{eqnarray*}
\frac{E(\De x)}{\De x} \stackrel{\De x\ra 0}{\ra}\quad\mbox{const}.
\end{eqnarray*}
$$(ref4)

Now the difference quotient {eq}`ForwardDiff` is not the only
one that can be used to approximate $f'(x)$. Indeed, if we consider
the expansion of $f(x-\De x)$:

$$
\begin{eqnarray}
f(x-\De x)=f(x)-\De xf'(x)+\frac{\De x^2}{2!}f''(x)-\frac{\De
x^3}{3!}f^{(3)}( x)+\frac{\De x^4}{4!}f^{(4)}(x)+\ldots
\label{Taylorminus}
\end{eqnarray}
$$(Taylorminus)

and we subtract {eq}`Taylorminus` from {eq}`Taylor`, and divide
by $(2\De x)$ we obtain:

$$
\begin{eqnarray}
\frac{f(x+\De x)-f(x-\De x)}{2\De x}=f'(x)+\frac{\De
x^2}{3!}f^{(3)}(\xi ) \label{CentralDiff1}
\end{eqnarray}
$$(CentralDiff1)

We notice that the error term associated with this form of
difference approximation is $O(\De x^2)$, which converges more
rapidly to zero as $\De x\ra 0$.

In order to obtain an approximation to $f''(x)$ we add
{eq}`Taylorminus` to {eq}`Taylor`, which upon re-arrangement and
dividing by $\De x^2$ leads to:

$$
\begin{eqnarray}
\frac{f(x+\De x)-2f(x)+f(x-\De x)}{\De x^2}=f''(x)+\frac{1}{12}\De
x^2f^{(4)}(\xi ) \label{CentralDiff2}
\end{eqnarray}
$$(CentralDiff2)

Due to the symmetry of the difference approximations
{eq}`CentralDiff1` and {eq}`CentralDiff2` about the expansion
point $x$ these are called central difference approximations. The
difference approximation {eq}`ForwardDiff` is known as a forward
difference approximation. We note that the central difference
schemes {eq}`CentralDiff1` and {eq}`CentralDiff2` are second
order accurate while the forward difference scheme
{eq}`ForwardDiff` is only accurate to $O(\De x)$.

### Solving The Heat Equations Using The Method Of Finite Differences

Consider the following initial-boundary value problem for the heat
equation

$$
\begin{eqnarray}
\frac{\pd u}{\pd t} &=& \al^2\frac{\pd ^2u}{\pd x^2}\quad 0<x<1, t>0 \label{HeatEq}\\
\mbox{BC:}\quad u(0,t)&=& 0\quad u(1,t)=0 \label{DirichletBC}\\
\mbox{IC:}\quad u(x,0) &=& f(x).
\end{eqnarray}
$$(HeatEq)

The basic idea is to replace the derivatives in the heat equation by
difference quotients. We consider the relationships between $u$ at
$(x,t)$ and its neighbours a distance $\De x$ apart and at a time
$\De t$ later.

Corresponding to the difference quotient approximations introduced
in {ref}`fd-intro`, we consider the following partial difference
approximations.

___Forward Difference in Time:___

$$
\begin{eqnarray*}
u(x,t+\De t)=u(x,t)+\De t\frac{\pd u}{\pd t} (x,t)+\frac{\De
t^2}{2!}\frac{\pd ^2u}{\pd t^2} (x,t)+\cdots .
\end{eqnarray*}
$$(ref9)

After re-arrangement and division by $\De t$:

$$
\begin{eqnarray}
\frac{u(x,t+\De t)-u(x,t)}{\De t} =\frac{\pd u}{\pd t}(u,t)+O(\De t)
\label{TimeDifference}
\end{eqnarray}
$$(TimeDifference)

___Central Differences in Space:___

$$
\begin{eqnarray*}
u(x+\De x,t) &=& u(x,t)+\De x\frac{\pd u}{\pd x}(x,t) +\frac{\De
x^2}{2!}\frac{\pd^2u}{\pd x^2}(u,t) +\frac{\De x^3}{3!}\frac{\pd
^3u}{\pd x^3}(x,t)
+\frac{\De x^4}{4!}\frac{\pd ^4u}{\pd x^2}(x,t)+\cdots \\
u(x-\De x,t) &=& u(x,t)-\De x\frac{\pd u}{\pd x}(x,t) +\frac{\De
x^2}{2!}\frac{\pd ^2u}{\pd x^2}(x,t) -\frac{\De x^3}{3!}\frac{\pd
^3u}{\pd x^3}(x,t) +\frac{\De x^4}{4!}\frac{\pd ^4u}{\pd
x^4}(x,t)+\cdots .
\end{eqnarray*}
$$(ref11)

Adding and re-arranging:

$$
\begin{eqnarray}
\frac{u(x+\De x,t)-2u(x,t)+u(x-\De x,t)}{\De x^2}=\frac{\pd ^2u}{\pd
x^2}(x,t) +O(\De x^2) \label{CentralDiffSpace}
\end{eqnarray}
$$(CentralDiffSpace)

Substituting {eq}`TimeDifference` and {eq}`CentralDiffSpace`
into {eq}`HeatEq` we obtain

$$
\begin{eqnarray*}
\frac{u(x,t+\De t)-u(x,t)}{\De t}=\al ^2\left( \frac{u(x+\De
x,t)-2u(x,t)+u(x-\De x,t)}{\De x^2}\right) +O(\De t,\De x^2) .
\end{eqnarray*}
$$(ref13)

Re-arranging:

$$
\begin{eqnarray}
u(x,t+\De t)=u(x,t)+\al ^2\left(\frac{\De t}{\De x^2}\right) \left\{
u(x+\De x,t)-2u(x,t)+u(x-\De x,t)\right\} .
\end{eqnarray}
$$(ref14)

We subdivide the spatial interval $[0,1]$ into $N+1$ equally spaced
sample points $x_n=n\De x$. The time interval $[0,T]$ is subdivided
into $M+1$ equal time levels $t_k=k\De t$. At each of these
space-time sample points we introduce approximations:

$$
\begin{eqnarray*}
u(x_n,t_k)\simeq u_n^k.
\end{eqnarray*}
$$(ref15)

```{figure} ../img/fd/heat_eq_scheme.png
:name: fd_scheme
:alt: fd_schem
:align: center
```

$$
\begin{eqnarray*}
u_n^{k+1}=u_n^k+\al ^2\left(\frac{\De t}{\De x^2}\right) \left(
u_{n+1}^k-2u_n^k+u_{n-1}^k\right)
\end{eqnarray*}
$$(ref16)

___Implementing Derivative Boundary Conditions:___

Assume that the boundary conditions for {eq}`HeatEq` are changed
to

$$
\begin{eqnarray*}
\mbox{BC:}\quad u(0,t)=0,\quad\frac{\pd u}{\pd x}(1,t)=0.
\end{eqnarray*}
$$(ref17)

Consider a central difference approximation to $\dis\frac{\pd u}{\pd
x}(1,t)$, where $x_N=N \Delta x = 1$,

$$
\begin{eqnarray*}
\frac{u(x_N+\De x,t)-u(x_N-\De x,t)}{\De x}=0.
\end{eqnarray*}
$$(ref18)

Re-arranging we obtain:

$$
\begin{eqnarray*}
u(x_N+\De x,t)=u(x_N-\De x,t)\quad (*)
\end{eqnarray*}
$$(ref19)

Since $x_N=1$ we observe that $x_N+\De x$ is outside the domain. To
accomodate this we introduce an extra column $u_{N+1}$ into which we
copy the values $u_{N-1}$. In the column $x_N$ we implement the same
difference approximation for the Heat Equation, namely:

$$
\begin{eqnarray*}
u_N^{k+1}=u_N^k+\al ^2( \frac{\De t}{\De x^2} )
(u_{N+1}^k-2u_N^k+u_{N-1}^k) \quad (**)
\end{eqnarray*}
$$(ref20)

While $u_{N+1}^k=u_{N-1}^k $ (see (*) ) since column
$u_{N-1}^k$ is copied to column $u_{N+1}^k$. 

```{note}
This BC could be implemented another way without introducing the additional
column, by eliminating $u_{N+1}$ from $(*)$ and $(**)$:

$$
\begin{eqnarray*}
u_{N}^{k+1}=u_N^k+2\al ^2\left(\frac{\De t}{\De x^2}\right)\left(
u_{N-1}^k-u_N^k\right) .
\end{eqnarray*}
$$(ref21)

If this latter equation is implemented at $x_N$ there is no need to
introduce an extra column $u_{N+1}$ or to implement the difference
equation given in (**) as the the derivative boundary condition is
taken care of automatically.
```

#### Stability of the Finite Difference Scheme for the Heat Equation  

Consider the following finite difference approximation to the 1D heat equation:

$$
\begin{eqnarray*}
u_n^{k+1}-u_n^k=\frac{\Delta t}{\Delta x^2} \left(
u_{n+1}^k-2u_n^k+u_{n-1}^k\right)\quad \mbox{where} \quad u_n^k\simeq
u(x_n,t_k)
\end{eqnarray*}
$$(ref22)
Let $\dis u_n^k=\phi_ke^{in\De x\theta}$ then

$$
\begin{eqnarray*}
(\phi_{k+1}-\phi_k) e^{in\De x\theta} &=& \frac{\De t}{\Delta x^2} \left( e^{i\De x\theta}-2+e^{-i\De x\theta}\right) \phi_k e^{in\De x\theta}\\
&=& \frac{\De t}{\De x^2}\left[ 2\cos (\theta\De x)-2\right]\phi_k e^{in\De x\theta}
\end{eqnarray*}
$$(ref23)

Therefore

$$
\begin{eqnarray*}
\phi_{k+1}&=&\phi_k-\frac{\De t}{\De x^2}4\sin^2\left(\frac{\theta\De x}{2}\right)\phi_k \quad \mbox{since} \quad\cos (\theta\De x)-1=-2\sin^2 \left(\frac{\theta\De x}{2}\right)\\
&=&\left[1-\frac{4\De t}{\De x^2}\sin^2\left(\frac{\theta\De
x}{2}\right)\right]\phi_k
\end{eqnarray*}
$$(ref24)

Now for stability we require that $|\phi_{k+1}|\leq |\phi_k|$ so that

$$
\begin{eqnarray*}
&&\left| 1-\frac{4\De t}{\De x^2}\sin^2\left(\frac{\theta\De x}{2}\right)\right|\leq 1\\
\rightarrow &&-2\leq -\frac{4\De t}{\De
x^2}\sin^2\left(\frac{\theta\De x}{2}\right)\leq 0
\end{eqnarray*}
$$(ref25)

The right inequality is satisfied automatically, while the left inequality can be re-written 
in the form:

$$
\begin{eqnarray*}
&\frac{4\De t}{\De x^2}\sin^2\left(\frac{\theta\De x}{2}\right)\leq 2
\end{eqnarray*}
$$(ref26)

Since $\sin^2(\ \ )\leq 1$ this condition is satisfied for all
$\theta$ provided

$$
\begin{eqnarray}
\De t\le\frac{\De x^2}{2} \label{eqStabilityCondHeat}
\end{eqnarray}
$$(eqStabilityCondHeat)

#### Exercises on Finite Differences Applied to the Heat Equation

___Numerical Instability:___

 (a) Change the $\De t$ in cell D1 from $0{.}001$ to $0{.}05$ and
you will observe what is known as a numerical instability. Now
change $\De t$ to $0{.}00625$, which is known as the stability
boundary predicted by {eq}`eqStabilityCondHeat` and observe what
happens. Now let $\De t=0{.}006$ and observe the abrupt change in
the solution - it is much closer to what we would expect.

 (b) Derive the stability condition for
the finite difference approximation of the 1D heat equation when
$\alpha^2\ne 1$.

$$
\begin{eqnarray*}
\mbox{i.e.}\quad u_n^{k+1}-u_n^k=\frac{\alpha^2\De t}{\De x^2}\left(
u_{n+1}^k-2u_n^k+u_{n-1}^k\right)
\end{eqnarray*}
$$(ref28)

___Truncation Error:___

The instability noted in 1 above is not the only source of error in the numerical approximation.
Although numerical instability is evident for a parameter choice
that is unstable, the other type of error is present in almost every
type of numerical approximation scheme.  This class of error results
from discarding the $O(\De x^2)$ and $O(\De t)$ terms in {eq}`CentralDiff1` and {eq}`ForwardDiff`
when we replace derivatives in {eq}`HeatEq` by difference quotients. This error is known as the truncation error. To estimate the magnitude of the truncation
error, change the spread sheet to implement the initial condition

$$
\begin{eqnarray*}
f(x)=\left\{\begin{array}{ll}2x &0<x<1/2\\ 2(1-x)&1/2\leq
x<1\end{array}\right. .
\end{eqnarray*}
$$(ref29)

Now code up the Fourier Series (in another spread sheet) that is
derived in Lecture 10, Exercise 10.1 and compare the numerical
solution to the `exact' Fourier Series solution with 50 terms. The
difference between the two is mainly due to the truncation error
since the round-off error is about $10^{-12}$ and does not grow if
stable parameters are used.

___Derivative Boundary Conditions:___

Implement a derivative boundary condition the left endpoint $x=0$.
Check the numerical solution against the problem solved in Lecture
11.

### Finite Difference Scheme For The 1D Wave Equation

Consider the following initial boundary value problem for the Wave
Equation:

$$
\begin{eqnarray}
u_{tt}&=&c^2u_{xx}\quad 0<x<L\\
\mbox{BC:}\quad u(0,t) &=& 0\quad u(L,t)=0\\
\mbox{IC:}\quad u(x,0)&=&f(x)\\
\frac{\pd u}{\pd t}(x,0)&=&g(x) \label{eqwaveVelIC}
\end{eqnarray}
$$(eqwaveVelIC)

We introduce a finite difference mesh $x_n=n\De t$, $t_k=k\De t$ and
let the corresponding nodal values be denoted by

$$
\begin{eqnarray*}
u_n^k\simeq u(x_n,t_k).
\end{eqnarray*}
$$(ref31)

```{figure} ../img/fd/wave_eq_scheme.png
:name: fd_scheme_wave
:alt: fd_scheme_wave
:align: center
```

Now approximating derivatives by central differences both in space
and time we obtain

$$
\begin{eqnarray}
\frac{u_n^{k+1}-2u_n^k+u_n^{k-1}}{\De t^2}&=&c^2\left(\frac{u_{n+1}^k-2u_n^k+u_{n-1}^k}{\De x^2}\right) +O(\De x^2,\De t^2).\\
\mbox{Therefore}\quad u_n^{k+1}&=&2u_n^k-u_n^{k-1}+{\left(\frac{c\De t}{\De x}\right)}^2\big( u_{n+1}^k-2u_n^k+u_{n-1}^k\big)\nonumber\\
\underbrace{u_n^{k+1}}_{\mbox{time level $k+1$}}&=&
\underbrace{r^2u_{n+1}^k+2(1-r^2)u_n^k+r^2u_{n-1}^k}_{\mbox{time
level $k$}}-\underbrace{u_n^{k-1}}_{\mbox{time level $k-1$}}
\label{eq_FD_Wavekgen}
\end{eqnarray}
$$(eq_FD_Wavekgen)

Here $r=(c\De t/\De x)$ is known as the Courant Number. We observe
that the Discrete Equation {eq}`eq_FD_Wavekgen` involves three
distinct levels in which known data is transferred from steps $k-1$
and $k$ to step $k+1$.

___Initial Conditions - Starting the Solution___

The 3-level scheme poses some challenges when imposing the initial
conditions. If we imagine a row of false mesh points at time $t=-\De
t=t_{-1}$, then the initial velocity condition {eq}`eqwaveVelIC`
can be approximated using central differences as:

$$
\begin{eqnarray}
\frac{u_n^1-u_n^{-1}}{2\De t}& = & g(x_n)\label{FalseICWaveEq}
\end{eqnarray}
$$(ref33)

therefore

$$
\begin{eqnarray}
u_n^{-1}&=&u_n^1-2 \Delta t g(x_n) \label{eqWaveIC}
\end{eqnarray}
$$(eqWaveIC)

Now we assume that the Discrete Wave Equation {eq}`eq_FD_Wavekgen`
also holds at $t=0$ so that

$$
\begin{eqnarray}
u_n^1=r^2u_{n+1}^0+2(1-r^2)u_n^0+r^2u_{n-1}^0-u_n^{-1}
\label{eq_FD_Wavekeq1}
\end{eqnarray}
$$(eq_FD_Wavekeq1)

Substituting {eq}`eqWaveIC` into {eq}`eq_FD_Wavekeq1` and
re-arranging we obtain:

$$
\begin{eqnarray}
u_n^1=\frac{1}{2}(r^2u_{n+1}^0+2(1-r^2)u_n^0+r^2u_{n-1}^0)+\De t
g(x_n) \label{eq_FD_Wavekeq1c}
\end{eqnarray}
$$(eq_FD_Wavekeq1c)

Since $u_n^0=f(x_n)$ and $g(x_n)$ are known, we are now in a
position to specify the first two rows of nodes. This is sufficient
to start the recursion {eq}`eq_FD_Wavekgen` for all subsequent
steps.

```{figure} ../img/fd/wave_scheme_false_mesh.png
:name: false_mesh
:alt: wave_scheme_false_mesh
:align: center

8.20 refers to {eq}`eq_FD_Wavekgen` and 8.24 refers to {eq}`eq_FD_Wavekeq1c`.
$\bullet=$ nodes and $\circ=$ false mesh points used to derive {eq}`eq_FD_Wavekeq1c`
but not actually used in the computation.
```

#### Stability of the Finite Difference Scheme for the Wave Equation

Consider the following finite difference approximation to the 1D
wave equation:

$$
\begin{eqnarray}
u_n^{k+1}=r^2u_{n+1}^k+2(1-r^2)u_n^k+r^2u_{n-1}^k-u_n^{k-1}
\label{FDWaveEqGen}
\end{eqnarray}
$$(FDWaveEqGen)

and, as in the case of the heat equation, substitute $\dis
u_n^k=\phi_ke^{in\De x\theta}$ into {eq}`FDWaveEqGen`

$$
\begin{eqnarray*}
e^{in\Delta x\theta }\phi _{k+1}=\left( r^{2}e^{i\Delta x\theta
}+2\left( 1-r^{2}\right) +r^{2}e^{-i\Delta x\theta }\right)
e^{in\Delta x\theta }\phi _{k}-e^{in\Delta x\theta }\phi _{k-1}
\end{eqnarray*}
$$(ref38)

Canceling terms and using the double angle formulae

$$
\begin{eqnarray}
\phi _{k+1} &=&2\left( 1+r^{2}\left( \cos \Delta x\theta -1\right)
\right)
\phi _{k}-\phi _{k-1} \\
&=&2\left( 1-2r^{2}\sin ^{2}\frac{\Delta x\theta }{2}\right) \phi
_{k}-\phi _{k-1} \label{waveStabEq}
\end{eqnarray}
$$(waveStabEq)

 If we now assume that $\phi _{k}$ has the following exponential
form $\phi _{k}=G^{k}$ then {eq}`waveStabEq` reduces to the
following quadratic equation:

$$
\begin{eqnarray}
G^{2}-2\gamma G+1=0 \label{StabQuadratic}
\end{eqnarray}
$$(ref40)

where $\gamma =\left( 1-2r^{2}\sin ^{2}\frac{\Delta x\theta
}{2}\right) $. The solutions of this quadratic equation are given by

$$
\begin{eqnarray}
G_{1,2}=\gamma \pm \sqrt{\gamma ^{2}-1}
\end{eqnarray}
$$(ref41)

Now since $G_1$ and $G_2$ are the roots of this quadratic we may
conclude that

$$
\begin{eqnarray}
(G-G_{1})(G-G_{2})=G^{2}-(G_{1}+G_{2})G+G_{1}G_{2}=0
\label{StabQuadratic1}
\end{eqnarray}
$$(StabQuadratic1)

Comparing the last terms in these two quadratic equations
(\ref{StabQuadratic}) and {eq}`StabQuadratic1` we conclude that

$$
\begin{eqnarray}
G_{1}G_{2}=1. \label{unitDiskConstraint}
\end{eqnarray}
$$(unitDiskConstraint)

However, for stability of solutions for the form $\phi _{k}=G^{k}$,
we require that $|G_1|\le 1$ and $|G_2|\le 1$. Given the constraint
{eq}`unitDiskConstraint`, the only possibility, if the solutions
are to be stable, is that $|G_1|=|G_2|= 1$. Thus $G$ must fall on
the unit disk, which implies that

$$
\begin{eqnarray*}
|\gamma|\le 1
\end{eqnarray*}
$$(ref44)

Thus

$$
\begin{eqnarray*}
\left| 1-2r^{2}\sin ^{2}\frac{\Delta x\theta }{2}\right|\le 1
\end{eqnarray*}
$$(ref45)

or

$$
\begin{eqnarray*}
-1\le  1-2r^{2}\sin ^{2}\frac{\Delta x\theta }{2}\le 1
\end{eqnarray*}
$$(ref46)

so that

$$
\begin{eqnarray}
-2\le -2r^{2}\sin ^{2}\frac{\Delta x\theta }{2}\le 0
\label{StabCondWaveEq}
\end{eqnarray}
$$(StabCondWaveEq)

The second inequality in {eq}`StabCondWaveEq` is satisfied
automatically, while the first leads to the condition

$$
\begin{eqnarray*}
r^{2}\sin ^{2}\frac{\Delta x\theta }{2}\le 1
\end{eqnarray*}
$$(ref48)

Since the maximum value that $\sin ^{2}(\frac{\Delta x\theta }{2})$
can achieve is 1, we conclude that the condition for stability is

$$
\begin{eqnarray*}
r=(c\De t/\De x)\le 1
\end{eqnarray*}
$$(ref49)

or

$$
\begin{eqnarray}
\De t\le  \frac{\Delta x }{c} \label{CFL_Wave}
\end{eqnarray}
$$(CFL_Wave)

The condition {eq}`CFL_Wave`, which imposes an upper bound on the
time step that can be used, is known as the Courant-Friedrichs-Lewy
or CFL condition.

### Solving Laplace's Equation Using Finite Differences

Consider the boundary value problem

$$
\begin{eqnarray}
& & \frac{\dis\pd^2u}{\dis \pd x^2} +\frac{\pd ^2u}{\pd y^2}=0\quad 0<x,y<1 \label{Laplace}\\
\mbox{BC:}\quad & u(0,y)=0; &\quad u(1,y)=0;\quad u(x,0)=f(x);\quad
u(x,1)=0. \label{LaplaceBC}
\end{eqnarray}
$$(Laplace)

```{figure} ../img/fd/laplace_scheme.png
:name: laplace_scheme
:alt: laplace_scheme
:align: center
```

As before we replace the second derivatives in {eq}`Laplace` by
central difference quotients that are second order accurate:

$$
\begin{eqnarray}
\frac{u(x+\De x,y)-2u(x,y)+u(x-\De x,y)}{\De x^2} & = & \frac{\pd^2u}{\pd x^2}(x,y)+O(\De x^2) \label{CentDiffX}\\
\frac{u(x,y+\De y)-2u(x,y)+u(x,y-\De y)}{\De y^2} & = &
\frac{\pd^2u}{\pd y^2}(x,y)+O(\De y^2) \label{CentDiffY}
\end{eqnarray}
$$(CentDiffY)

We partition the interval $0\leq x\leq 1$ into $(N+1)$ equally
spaced nodes $x_n=n\De x$ and the interval $0\leq y\leq 1$ into
$(M+1)$ equally spaced nodes $y_m=m\De y$. Replacing the derivatives
in {eq}`Laplace` by the difference quotients in (\ref{CentDiffX})
and {eq}`CentDiffY`, and representing the mesh values at
$(x_n,y_m)$ by $u_{nm}\simeq u(x_n,y_m)$ we obtain:

$$
\begin{eqnarray*}
\frac{u_{n+1m}-2u_{nm}+u_{n-1m}}{\De
x^2}+\frac{u_{nm+1}-2u_{nm}+u_{nm-1}}{\De
y^2}={(u_{xx}+u_{yy})}_{(x_n,x_m)}+O(\De x^2,\De y^2).
\end{eqnarray*}
$$(ref53)

If we choose $\De x=\De y$ then we obtain

$$
\begin{eqnarray}
u_{n+1m}+u_{n-1m}+u_{nm+1}+u_{nm-1}-4u_{nm}=0\quad 1\leq n,m\leq
(N-1),(M-1) \label{DiscreteLaplace}
\end{eqnarray}
$$(DiscreteLaplace)

```{figure} ../img/fd/laplace_stencil.png
:name: laplace_stencil
:alt: laplace_stencil
:align: center
```

This is known as the finite difference `Stencil' that relates
$u_{nm}$ to its 4 nearest neighbours.

 This is a system of $(N-1)\times (M-1)$ unknowns for the values
of $u_{nm}$ interior to the domain - recall the boundary values are
already specified!

#### Solving the System of Equations by Jacobi Iteration

This is a procedure to solve the system of Equations
{eq}`DiscreteLaplace` by looping through each of the mesh points
and updating $u_{nm}$ according to {eq}`DiscreteLaplace` assuming
that the nearest neighbours already have values close to the exact
solution. This procedure is repeated until the changes that are made
in each iteration falls below a certain tolerance.

To implement this iterative procedure we observe that the discrete
Laplace Equation {eq}`DiscreteLaplace` can be re-written in the
form:

$$
\begin{eqnarray}
u_{nm}^{k+1}=\frac{u_{n+1m}^k+u_{n-1m}^k+u_{nm+1}^k+u_{nm-1}^k}{4}
\label{Jacobi}
\end{eqnarray}
$$(Jacobi)

```{figure} ../img/fd/laplace_average.png
:name: laplace_average
:alt: laplace_average
:align: center
```

Thus $u_{nm}$ is the average value of its nearest neighbours. Note
that a new superscript index $k$ has been introduced to represent
the nodal values at the $k$th iteration. Thus iteration can be
viewed as taking successive neighbour averages until there is no
change, at which point the value of $u_{mn}$ equals the average of
the values at its mesh neighbours. This mean value property is a
discrete form of a fundamental property of any solution to Laplace's
Equation, which also implies that maxima or minima of solutions to
Laplace's equation cannot occur inside the domain but must be
restricted to the boundary - a very useful poperty of Laplace's
equation known as The Maximum Principle.

To implement the iterative procedure {eq}`Jacobi` in a spread
sheet, go to the Tools Menu at the top of the screen and click
on the Options Tab. Then select the Calculation Tab.
Check the Iteration box. If you set the number of iterations
to 5 say, then if you start with zero values throughout the interior
of the domain (as you should if you cut and paste as demonstrated in
class), you will see the values percolate 5 cells into the domain
from the non zero boundary condition $f(x)=\sin (\pi x)$. You can
choose a surface plot to visualize the solution. Now hold down the
F9 key and watch the solution move to equilibrium. This
iterative process essentially uses diffusion on a pseudo time scale
to take the solution to equilibrium.\medskip

#### Exercises for Laplace's Equation

1. Implement a 0 derivative BC along the lines $x=0$ and $x=1$:

Plot a cross section of the results along $y=1/2$. To ensure that
$\dis\frac{\pd u}{\pd x}(0,y)=0=\frac{\pd u}{\pd x}(1,y)$

2. Implement an inhomogeneous term for Poisson's Equation:

$$
\begin{eqnarray*}
\frac{\pd^2u}{\pd x^2}+\frac{\pd^2u}{\pd y^2}=f(x,y)\quad 0<x,y<1.
\end{eqnarray*}
$$(ref56)

Introduce finite difference quotients, assume $\De x=\De y$ to
arrive at the iterative formula:

$$
\begin{eqnarray*}
u_{nm}^{k+1}=\frac{\big(u_{n+1m}^k+u_{n-1m}^k+u_{nm+1}^k+u_{nm-1}^k-\De
x^2f(x_n,y_m)\big)}{4}.\quad (*)
\end{eqnarray*}
$$(ref57)

It may be useful to calculate the values of $f_{nm}$ on a separate
sheet in which the same cell values as those for $u_{nm}$ are
maintained. Then the values of $f_{nm}$ can be referenced in the
calculation of $u_{nm}$ according to $(*)$.
