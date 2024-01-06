# Laplace Equation

Consider the boundary value problem

$$
\begin{eqnarray}
& & \frac{\displaystyle\partial^2u}{\displaystyle \partial x^2} +\frac{\partial ^2u}{\partial y^2}=0\quad 0<x,y<1 \label{Laplace}\\
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
\frac{u(x+\Delta x,y)-2u(x,y)+u(x-\Delta x,y)}{\Delta x^2} & = & \frac{\partial^2u}{\partial x^2}(x,y)+O(\Delta x^2) \label{CentDiffX}\\
\frac{u(x,y+\Delta y)-2u(x,y)+u(x,y-\Delta y)}{\Delta y^2} & = &
\frac{\partial^2u}{\partial y^2}(x,y)+O(\Delta y^2) \label{CentDiffY}
\end{eqnarray}
$$(CentDiffY)

We partition the interval $0\leq x\leq 1$ into $(N+1)$ equally
spaced nodes $x_n=n\Delta x$ and the interval $0\leq y\leq 1$ into
$(M+1)$ equally spaced nodes $y_m=m\Delta y$. Replacing the derivatives
in {eq}`Laplace` by the difference quotients in (\ref{CentDiffX})
and {eq}`CentDiffY`, and representing the mesh values at
$(x_n,y_m)$ by $u_{nm}\simeq u(x_n,y_m)$ we obtain:

$$
\begin{eqnarray*}
\frac{u_{n+1m}-2u_{nm}+u_{n-1m}}{\Delta
x^2}+\frac{u_{nm+1}-2u_{nm}+u_{nm-1}}{\Delta
y^2}={(u_{xx}+u_{yy})}_{(x_n,x_m)}+O(\Delta x^2,\Delta y^2).
\end{eqnarray*}
$$(ref-fd-pdes-53)

If we choose $\Delta x=\Delta y$ then we obtain

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
$\displaystyle\frac{\partial u}{\partial x}(0,y)=0=\frac{\partial u}{\partial x}(1,y)$

2. Implement an inhomogeneous term for Poisson's Equation:

$$
\begin{eqnarray*}
\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}=f(x,y)\quad 0<x,y<1.
\end{eqnarray*}
$$(ref-fd-pdes-56)

Introduce finite difference quotients, assume $\Delta x=\Delta y$ to
arrive at the iterative formula:

$$
\begin{eqnarray*}
u_{nm}^{k+1}=\frac{\big(u_{n+1m}^k+u_{n-1m}^k+u_{nm+1}^k+u_{nm-1}^k-\Delta
x^2f(x_n,y_m)\big)}{4}.\quad (*)
\end{eqnarray*}
$$(ref-fd-pdes-57)

It may be useful to calculate the values of $f_{nm}$ on a separate
sheet in which the same cell values as those for $u_{nm}$ are
maintained. Then the values of $f_{nm}$ can be referenced in the
calculation of $u_{nm}$ according to $(*)$.
