# Separation of Variables and Fourier Series

In this lecture we will introduce the method of separation of variables by using it to solve the heat equation, which reduces the solution of the PDE to solving two ODEs, one in time and one in space. The time ODE represents the exponential decay of the solution, while the second order spatial ODE along with the boundary conditions give rise to an eigenvalue problem, which needs to be solved to identify the so-called "separation" constant. Super-position of the resulting solutions leads naturally to the expansion of the initial temperature distribution $f(x)$ in terms of a series of sine functions - known as a Fourier Series.

```{admonition} Key Concepts
Heat equation; boundary conditions; Separation of variables; Eigenvalue problems for ODE; Fourier Series.
```

## The Heat/Diffusion equation and dispersion relation

We consider the heat equation (or diffusion equation)

$$
\frac{\partial u}{\partial t} = \alpha^2 \frac{\partial^2 u}{\partial x^2} 
$$ (eq:heat)

where $\alpha^2$ is the thermal conductivity. If we look for exponential solutions of the form

$$
u(x,t) = e^{i k x + \sigma t}
$$ (eq:ExponentialSpaceTimeSolution)

we obtain the dispersion relation $\sigma=-k^2 \alpha^2$ and corresponding solutions

$$
u(x,t) = e^{-k^2 \alpha^2 t} e^{i k x}
$$ (eq:ExponentialSpaceTimeSolutionafterdispersionRel)

We observe that in this family of solutions, which are parameterized by $k$, the solution is decomposed into the product of a time function that decays exponentially with increasing $t$ and a spatial function comprising $\sin(kx)$ and $\cos(kx)$. The parameter $k$, which is called the *wavenumber*, needs to be determined in
order to match the boundary conditions in the problem. We will see that for problem defined on a finite domain there are a countable infinity of admissible $k$ values, which we will determine by solving an eigenvalue problem involving an ordinary differential operator.

## Types of Initial-Boundary Value Problems for the heat equation

We consider the heat equation subject to the following initial and boundary conditions:

### Dirichlet Boundary Conditions

```{figure} ../img/fourier/DirichletBar.png
---
width: 80%
name: DirichletBar
---
Consider a conducting bar with thermal conductivity $\alpha^2$ that has an initial temperature distribution $u(x,0)=f(x)$ and whose endpoints are maintained at $0^{\circ }C$, i.e. embedded in ice
```

What do you expect the solution to look like as $t \to \infty$?

### Neumann Boundary Conditions

```{figure} ../img/fourier/NeumannBar.png
---
width: 80%
name: NeumannBar
---
Consider a conducting bar with thermal conductivity $\alpha^2$ that has an initial temperature distribution $u(x,0)=f(x)$ and whose endpoints are insulated
```

What do you expect the solution to look like as $t \to \infty$?

### Mixed Boundary Conditions

```{figure} ../img/fourier/MixedBar.png
---
width: 80%
name: MixedBar
---
Consider a conducting bar with thermal conductivity $\alpha^2$ that has an initial temperature distribution $u(x,0)=f(x)$ and whose left endpoint is held at $^{\circ }C$ (i.e., embedded in a block of ice) while the right endpoint is insulated
```

What do you expect the solution to look like as $t \to \infty$?

## Separation of Variables - Fourier sine Series

Consider the heat conduction in an insulated rod whose endpoints are held at zero degrees for all time and within which the initial temperature is given by $f(x)$ as shown in {numref}`Figure {number} <DirichletBar>`.

*Fourier's Guess:*

$$
\begin{align}
u(x,t) &= X(x)T(t) \\
u_t &= X(x) \dot{T}(t) = \alpha^2 u_{xx} = \alpha^2 X''(x) T(t)
\end{align}
$$ (eq:Separation1)

$\div \alpha^2 XT$:

$$
\frac{X''(x)}{X(x)} = \frac{\dot{T}(t)}{\alpha^2 T(t)} = \mbox{ Constant } = -\lambda^2.
$$ (eq:Separation2)

$->$

$$
\begin{align}
T'(t) &= -\alpha^2 \lambda^2 T(t) \\
\frac{dT}{T} &= -\alpha^2 \lambda^2 dt \\
\ln|T| &= -\alpha^2 \lambda^2t + c \\
T(t) &= D e^{-\alpha^2 \lambda^2 t}
\end{align}
$$ (eq:Separation3)

$x>$

$$
\begin{align}
X''(x) + \lambda^2 X(x) &= 0 \\
\Rightarrow X(x) &= c_1 e^{i \lambda x} + c_{2} e^{-i \lambda x} \\
&= A\sin \lambda x + B\cos \lambda x
\end{align}
$$ (eq:Separation4)

*Impose the boundary conditions:*

$$
\begin{array}{lcrclcl}
0&=&u(0,t)&=&X(0)T(t)&=&BT(t)\Rightarrow B=0\\
0&=&u(L,t)&=&X(L)T(t)&=&(A\sin\lambda L)T(t)
\end{array}
$$ (eq:Separation5)

Now we do not want the trivial solution so $A\not= 0$. Thus we look for values of $\lambda$ such that

$$
\sin\lambda L = 0 \Rightarrow \lambda = \left( \frac{n\pi}{L} \right) \quad n=1,2,\ldots
$$ (eq:Separation6)

Thus

$$
u_n(x,t) = e^{-\alpha^2 {\left( \frac{n\pi}{L} \right)}^2 t} \sin\left(\frac{n\pi x}{L} \right) \quad n=1,2,\ldots
$$ (eq:Separation7)

are all solutions of the heat equation $u_t = \alpha^2u_{xx}$. Since the heat equation is linear, a linear combination of solutions is again a solution. Thus the most general solution is

$$
u(x,t) = \sum_{n=1}^\infty b_n \sin\left(\frac{n \pi x}{L}\right) e^{\alpha^2{\left(\frac{n\pi}{L}\right)}^2 t}
$$ (eq:heatsolution)

What about the initial condition $u(x,0)=f(x)$? Given $f(x)$ we need to find the $b_n$ such that the infinite series of functions

$$
\sum_{n=1}^\infty b_n\sin\left(\frac{n\pi x}{L}\right)
$$ (eq:heatsolutioninitial)

agrees with $f$ on $[0,L]$. Note that $f(x)$ may not be periodic $f(x+2L) \not= f(x)$ but the series is periodic since

$$
\sin\left(\frac{n\pi}{L}\right)(x + 2L) = \sin\left(\frac{n\pi x}{L}\right)
$$

In fact they do agree on $[0,L]$ and are different elsewhere.