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

# Separation Of Variables And Fourier Series

In this lecture we will introduce the method of separation of
variables by using it to solve the heat equation, which reduces the
solution of the PDE to solving two ODEs, one in time and one in
space. The time ODE represents the exponential decay of the
solution, while the second order spatial ODE along with the boundary
conditions give rise to an eigenvalue problem, which needs to be
solved to identify the so-called "separation" constant.
Super-position of the resulting solutions leads naturally to the
expansion of the initial temperature distribution $f(x)$ in terms of
a series of $\sin$ functions - known as a Fourier Series.

```{admonition} Key Concepts
Heat equation; boundary conditions;
Separation of variables; Eigenvalue problems for ODE; Fourier
Series.
```

## The Heat/Diffusion Equation And Dispersion Relation

We consider the heat equation (or diffusion equation)

$$
\begin{equation}
\frac{\partial u}{\partial t}=\alpha^2 \frac{\partial^2 u}{\partial
x^2} \label{eq:heat}
\end{equation}
$$(ref-intro-separation-0)

where $\alpha^2$ is the thermal conductivity. If we look for
exponential solutions of the form

$$
\begin{equation} u(x,t)
= \{\rm\ e\}^{\{\rm\ i\} kx + \sigma t},
 \label{eq:ExponentialSpaceTimeSolution}
\end{equation}
$$(ref-intro-separation-1)

we obtain the dispersion relation $\sigma=-k^2 \alpha^2$ and corresponding solutions

$$
\begin{equation} u(x,t)
= \{\rm\ e\}^{- k^2 \alpha^2 t}\{\rm\ e\}^{\{\rm\ i\} k x }.
 \label{eq:ExponentialSpaceTimeSolutionafterdispersionRel}
\end{equation}
$$(ref-intro-separation-2)

We observe that in this family of solutions, which are
parameterized by $k$, the solution is decomposed into the product of
a time function that decays exponentially with increasing $t$ and a
spatial function comprising $\sin(kx)$ and $\cos(kx)$. The parameter
$k$, which is called the wavenumber, needs to be determined in
order to match the boundary conditions in the problem. We will see
that for problem defined on a finite domain there are a countable
infinity of admissible $k$ values, which we will determine by
solving an eigenvalue problem involving an ordinary differential
operator.

## Types of Initial-Boundary Value Problems for the Heat Equation

We consider the heat equation subject to the
following initial and boundary conditions:

### Dirichlet Boundary Conditions

```{figure} ../img/intro/DirichletBarBCs.png
:name: DirichletBarBCs
:alt: DirichletBarBCs
:align: center

Consider a conducting bar with thermal conductivity
$\alpha^2$ that has an initial temperature distribution
$u(x,0)=f(x)$ and whose endpoints are maintained at $^{\circ }C$,
i.e. embedded in ice.
```

What do you expect the solution to look like as $t\rightarrow\infty$?

### Neumann Boundary Conditions

```{figure} ../img/intro/NeumannBarBCs.png
:name: NeumannBarBCs
:alt: NeumannBarBCs
:align: center

Consider a conducting bar with thermal conductivity
$\alpha^2$ that has an initial temperature distribution
$u(x,0)=f(x)$ and whose endpoints are insulated.
```

What do you expect the solution to look like as $t\rightarrow\infty$?

### Mixed Boundary Conditions

```{figure} ../img/intro/MixedBarBCs.png
:name: MixedBarBCs
:alt: MixedBarBCs
:align: center

Consider a conducting bar with thermal conductivity
$\alpha^2$ that has an initial temperature distribution
$u(x,0)=f(x)$ and whose left endpoint is held at $^{\circ }C$ (i.e.,
embedded in a block of ice) while the right endpoint is insulated.
```

What do you expect the solution to look like as $t\rightarrow\infty$?

## Separation Of Variables - Fourier Sine Series

Consider the heat conduction in an insulated rod whose endpoints are
held at zero degrees for all time and within which the initial
temperature is given by $f(x)$ as shown in figure
\ref{fig:DirichletBar}.

__Fourier's Guess:__

$$
\begin{eqnarray}
u(x,t) & = & X(x)T(t)\\
u_t & = & X(x)\dot{T}(t)=\alpha^2
u_{xx}=\alpha^2X^{\prime\prime}(x)T(t)\nonumber
\end{eqnarray}
$$(ref-intro-separation-3)

$\div\alpha^2 XT$:

$$
\begin{equation}
\frac{X^{\prime\prime}(x)}{X(x)}=\frac{\dot{T}(t)}{\alpha^2
T(t)}=\mbox{ Constant }=-\lambda^2.
\end{equation}
$$(ref-intro-separation-4)

$->$

$$
\begin{eqnarray}\begin{array}{c}
\dot{T}(t)=-\alpha^2\lambda^2 T(t)\quad\displaystyle\frac{dT}{T}=-\alpha^2\lambda^2\, dt\\
\ln|T| =-\alpha^2\lambda^2t+c\\
T(t)=D\{\rm\ e\}^{-\alpha^2\lambda^2t}.\end{array}
\end{eqnarray}
$$(ref-intro-separation-5)

$x>$

$$
\begin{eqnarray}\begin{array}{c}
X^{\prime\prime}(x)+\lambda ^2X(x)=0\\
\mbox{Guess}\quad X(x)=\{\rm\ e\}^{rx}\Rightarrow
(r^2+\lambda^2)\{\rm\ e\}^{rx}=0\quad r=\pm\lambda i\end{array}
\end{eqnarray}
$$(ref-intro-separation-6)

$$
\begin{eqnarray}\begin{array}{ccl}
X &= &c_1\{\rm\ e\}^{i\lambda x}+c_{2?}\{\rm\ e\}^{-i\lambda x}\\
  &= &A\sin\lambda x+B\cos\lambda x.\end{array}
\end{eqnarray}
$$(ref-intro-separation-7)

_Impose the boundary conditions:_

$$
\begin{eqnarray}\begin{array}{lcrclcl}
0&=&u(0,t)&=&X(0)T(t)&=&BT(t)\Rightarrow B=0\\
0&=&u(L,t)&=&X(L)T(t)&=&(A\sin\lambda L)T(t).\end{array}
\end{eqnarray}
$$(ref-intro-separation-8)

Now we do not want the trivial solution so $A\not= 0$. Thus we look
for values of $\lambda$ such that

$$
\begin{equation}
\sin\lambda L=0\Rightarrow\lambda =\left(\frac{n\pi}{L}\right)\quad
n=1,2,\ldots .
\end{equation}
$$(ref-intro-separation-9)

$$
\begin{eqnarray}
& &\mbox{Thus}\quad u_n(x,t)=\{\rm\ e\}^{-\alpha^2{\left(\frac{n\pi}{L}\right)}^2t} \sin\left(\frac{n\pi x}{L}\right)\quad n=1,2,\ldots\nonumber\\
& &\mbox{are all solutions of $u_t=\alpha^2u_{xx}$}.
\end{eqnarray}
$$(fourier_sin_solutions)

Since {eq}`fourier_sin_solutions` is linear, a linear combination of
solutions is again a solution. Thus the most general solution is

$$
\begin{equation}
u(x,t)=\sum\limits_{n=1}^\infty b_n\sin\left(\frac{n\pi
x}{L}\right)\{\rm\ e\}^{-\alpha^2{\left(\frac{n\pi}{L}\right)}^2}t.
\end{equation}
$$(ref-intro-separation-11)

What about the initial condition $u(x,0)=f(x)$.

$$
\begin{equation}
u(x,0)=f(x)=\sum\limits_{n=1}^\infty b_n\sin\left(\frac{n\pi
x}{L}\right).
\end{equation}
$$(ref-intro-separation-12)

Given $f(x)$ we need to find the $b_n$ such that the infinite series
of functions $\displaystyle\sum b_n\sin\left(\frac{n\pi x}{L}\right)$ agrees
with $f$ on $[0,L]$. \vspace{2in}

__Question:__ $f(x)$ may not be periodic $f(x+2L)\not= f(x)$
but the series is periodic since
$\displaystyle\sin\left(\frac{n\pi}{L}\right) (x+2L)=\sin\left(\frac{n\pi
x}{L}\right)$.

__Answer:__ In fact they do agree on $[0,L]$ and are different
elsewhere.
