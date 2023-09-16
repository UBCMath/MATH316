
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

# Fourier Sine Series

In the last lecture we reduced the problem of solving the
initial-boundary value problem for the heat distribution along a conducting
rod to solving two ODEs, one in space and one in time. The spatial
ODE and boundary conditions lead to an eigenvalue problem, which
identifies a discrete set of wavenumbers $\{\lambda_n\}$ and
corresponding eigenfunctions $\{X_n(x)\}$ that satisfy both homogeneous boundary
conditions and the spatial ODE. Because the heat equation is linear, the general
solution of the heat equation is obtained by superimposing the
product of the eigensolutions and the corresponding solution of
the time ODE to obtain an infinite series. All that remains is that we determine the
expansion coefficients $b_n$ for the terms of this series. These are
obtained by letting $t=0$ in the general solution and equating the
resulting series to the initial value function $f(x)$. This is
known as a Fourier Series. This lecture deals with the procedure to
determine the Fourier coefficients $b_n$. Our approach is motivated
by the process introduced in Linear Algebra for projecting a vector onto a set of basis vectors.

```{admonition} Key Concepts
Fourier Sine Series; Vector
Projection; functions as infinite dimensional vectors;
orthogonality; Fourier Coefficients.
```

Observe that we have a new type of eigenvalue problem in which we
seek a nontrivial solution to the following boundary value problem

$$
\begin{eqnarray}
LX&=&-X''=\lambda^2 X \quad \mbox{or}\quad X''+\lambda^2 X=0 \nonumber\\
X(0)&=&0=X(L). \label{eigenvalueBVP}
\end{eqnarray}
$$(ref-fourier-sine-0)

Just as in the case with matrices we obtain a sequence of
eigenvalues $\{\lambda_n\}$. However, because of the infinite
dimensional nature of this eigenvalue the problem there are an
infinite number of eigenvalues:

$$
\begin{equation}
\lambda _n=\left(\frac{n\pi}{L}\right)\quad n=1,2,\ldots
\end{equation}
$$(ref-fourier-sine-1)

and corresponding eigenfunctions $\{X_n(x)\}$

$$
\begin{eqnarray}\begin{array}{c}
X_n(x)=\sin\lambda _nx=\sin\left(\displaystyle\frac{n\pi x}{L}\right) \quad
\mbox{or}\quad X_n(x)\in\left\{
  \sin\left(\displaystyle\frac{\pi x}{L}\right) ,\sin\left(\displaystyle\frac{2\pi x}{L}\right) ,
   \sin\left(\displaystyle\frac{3\pi x}{L}\right) ,\ldots
\right\} .\end{array}
\end{eqnarray}
$$(ref-fourier-sine-2)

In order complete the solution of the heat equation we need
to determine the coefficients $b_n$ such that

$$
\begin{equation}
u(x,0)=f(x)=\sum\limits_{n=1}^\infty b_n\sin\left(\frac{n\pi
x}{L}\right) \label{eq:SineSeriesExpansionoff}.
\end{equation}
$$(ref-fourier-sine-3)

````{prf:observation}
:label: observation-fourier-sine-0
- For symmetric matrices the eigenvalues are real - for the BVP the eigenvalues $\{\lambda_n(x)\}$ are also real.
- For symmetric matrices the eigenvectors form a basis - the eigenfunctions $\{X_n(x)\}$ are linearly independent.
````

## Euler's Column: The Buckling Load for a Beam - Perhaps the Oldest Eigenvalue Problem

Eigenvalue problems also arise independently without necessarily coming
from a PDE problem. Consider a beam that is subjected to an axial
load  $P$ applied  to its endpoints. Our experience tells us that as
we increase $P$ a critical load $P_c$ is reached at which the beam
starts to buckle.

```{figure} ../img/fourier/euler_beam.png
:name: euler_beam
:align: center

Buckling of Euler's column.
```

The Bernoulli-Euler Law: A beam constructed from a
material that is known to deform in such a way that the curvature
$\kappa$ is proportional to the bending moment $M$.

In particular,

$$
\begin{eqnarray*}
\kappa=\frac{y''(x)}{[1+(y')^2]}=cM=\frac{M}{EI}
\end{eqnarray*}
$$(ref-fourier-sine-4)

where $E\ =$ Young's modulus and $I\ =$ the moment of inertia of the beam. If the deflection of the
beam is small $(y')^2<< |y'|<< 1$ then  we can make the approximation

$$
\begin{eqnarray*}
\framebox[6\totalheight]{$\displaystyle
y''=\frac{M}{EI}$}\quad\mbox{The Bernoulli-Euler Law}
\end{eqnarray*}
$$(ref-fourier-sine-5)

When subject to an axial load $P$ as shown in figure \ref{Euler_Beam_Figurec}, the bending moment on the buckling beam is given by

$$
\begin{eqnarray*}
M(x)&=&-Py(x)\\
y''&=&-\frac{M}{EI}=-\frac{P}{EI}y=-k^2y\quad k^2=\frac{P}{EI}
\end{eqnarray*}
$$(ref-fourier-sine-6)

Thus determining the magnitude of $P=EI k^2$ for which the beam will
first buckle is reduced to solving the following eigenvalue problem:

$$
\begin{eqnarray*}
\left.\begin{array}{l}y''+k^2y=0\\
y(0)=0=y(L)\end{array}\right\}\quad\mbox{Eigenvalue Problem}
\end{eqnarray*}
$$(ref-fourier-sine-7)

$$
\begin{eqnarray*}
y(x)&=&A\cos kx+B\sin kx\\
y(0)&=&A=0\quad y(L)=B\sin kL=0\Rightarrow k_n=\frac{n\pi}{L}\quad
n=1,2,\ldots .
\end{eqnarray*}
$$(ref-fourier-sine-8)

We have eigenvalues $\displaystyle k_n=\frac{n\pi}{L}$ and
eigenfunctions $\displaystyle y_n(x)=\sin\left(\frac{n\pi
x}{L}\right)$. But

$$
\begin{eqnarray*}
k_n^2=\frac{P_n}{EI}={\left(\frac{n\pi}{L}\right)}^2\quad
n=1,2,\ldots .
\end{eqnarray*}
$$(ref-fourier-sine-9)

Therefore the smallest buckling force $\displaystyle
P_c=P_1=\frac{EI\pi^2}{L^2}$ is known as the critical Euler load and
the Euler Buckling Mode is $\displaystyle\sin\left(\frac{\pi
x}{L}\right)$.

## Finding the Fourier Coefficients

How do we find the $b_n$ in the sine series expansion (\ref{eq:SineSeriesExpansionoff}) of $f(x)$?

```{figure} ../img/fourier/projection.png
:name: projection
:align: center

Left side: $\mbox{Expand } \mathbf{f}\mbox{ in terms of the
basis vectors}\big\{\mathbf{v}_1,\mathbf{v}_2,\mathbf{v}_3\big\}$.
Right side: Approximation of a function $f(x)$ by a vector
$\mathbf{f}$ of sample points with values $\{ f(x_k\}$.
```

**Decomposition of vectors into components - projection:**

How do we expand a vector $\mathbf{f}$ in terms of linearly
independent vectors $\mathbf{v}_k$?

$$
\begin{eqnarray}\begin{array}{l}
\mbox{Assume}\ \mathbf{f}= \al_1 \mathbf{v}_1 +\al_2 \mathbf{v}_2 +\al_3 \mathbf{v}_3\\
\mathbf{f}\cdot\mathbf{v}_k =\al_1 \mathbf{v}_1\cdot \mathbf{v}_k +\al_2 \mathbf{v}_2 \cdot\mathbf{v}_k +\al_3 \mathbf{v}_3\cdot \mathbf{v}_k\\
\quad\left[\begin{array}{ccc}
\mathbf{v}_1\cdot\mathbf{v}_1&\mathbf{v}_1\cdot\mathbf{v}_2&\mathbf{v}_1\cdot\mathbf{v}_3\\
\mathbf{v}_1\cdot\mathbf{v}_2          &\mathbf{v}_2\cdot\mathbf{v}_2          &\mathbf{v}_2\cdot\mathbf{v}_3\\
\mathbf{v}_1\cdot\mathbf{v}_3          &\mathbf{v}_2\cdot\mathbf{v}_3          &\mathbf{v}_3\cdot\mathbf{v}_3\\
\end{array}\right]\quad
\left[\begin{array}{c} \al_1 \\ \al_2 \\
\al_3\end{array}\right]\quad = \quad\left[\begin{array}{c}
\mathbf{f}\cdot\mathbf{v}_1\\
\mathbf{f}\cdot\mathbf{v}_2\\
\mathbf{f}\cdot\mathbf{v}_3\end{array}\right] \end{array}
\end{eqnarray}
$$(ref-fourier-sine-10)

If $\mathbf{v}_k\perp \mathbf{v}_\ell$, $k\not=\ell$ i.e. the
$\mathbf{v}_k$ are orthogonal

$$
\begin{equation}
\alpha _k=\frac{\mathbf{f}\cdot \mathbf{v}_k}{\mathbf{v}_k\cdot
\mathbf{v}_k}
\end{equation}
$$(ref-fourier-sine-11)

**Functions as infinite dimensional vectors and projection:**

But functions are just infinite dimensional vectors:

$$
\begin{eqnarray}
\mathbf{f} & \simeq & [f_1,f_2,\ldots ,f_N]\nonumber\\
\mathbf{g} & \simeq & [g_1,g_2,\ldots , g_N]\nonumber\\
\mathbf{f}\cdot \mathbf{g} & = & f_1g_1+f_2g_2+\cdots +f_Ng_N\qquad \Delta x=\frac{L}{N}\\
& = & \sum\limits_{k=1}^N f(x_k)g(x_k).\nonumber
\end{eqnarray}
$$(ref-fourier-sine-12)

The analogue of the dot product for functions is given by the
so-called inner product:

$$
\begin{equation}
\langle f,g\rangle :=\int\limits_0^L f(x)g(x)\, dx \simeq
\sum\limits_{k=1}^Nf(x_k)g(x_k)\Delta x.
\end{equation}
$$(ref-fourier-sine-13)

Back to finding $b_n$:

$$
\begin{eqnarray}
& & f(x)=\sum\limits_{n=1}^\infty b_n\sin\left(\frac{n\pi x}{L}\right)\\
& &\langle f,\sin\left(\frac{k\pi x}{L}\right)\rangle=
\int\limits_0^L f(x)\sin\left(\frac{k\pi x}{L}\right)\,
dx=\sum\limits_{n=1}^\infty b_n\int\limits_0^L \sin\left(\frac{n\pi
x}{L}\right)\sin\left(\frac{k\pi x}{L}\right)\, dx.\nonumber
\end{eqnarray}
$$(ref-fourier-sine-14)

Recall $\displaystyle\sin (A)\sin B=\frac{1}{2}\left\{\cos (A-B)-\cos
(A+B)\right\}$. Therefore

$$
\begin{eqnarray}
I_{nk} & = & \int\limits_0^L\sin\left(\frac{n\pi x}{L}\right)\sin\left(\frac{k\pi x}{L}\right)\, dx\nonumber\\
& = & \frac{1}{2}\int\limits_0^L\cos (n-k)\frac{\pi x}{L}-\cos (n+k)\frac{\pi x}{L}\, dx\qquad n\not= k\nonumber\\
& = & \frac{1}{2}{\left[\frac{\sin (n-k)\pi x/L}{(n-k)\pi /L}-\frac{\sin (n+k)\pi x/L}{(n+k)\pi /L}\right]}_0^L\nonumber\\
& = & 0\\
I_{nn} & = & \int\limits_0^L\sin^2 \left(\frac{n\pi x}{L}\right)\, dx=\frac{1}{2}\int\limits_0^L 1-\cos\left(\frac{2n\pi x}{L}\right)\, dx\nonumber\\
& = & L/2\nonumber
\end{eqnarray}
$$(ref-fourier-sine-15)

Therefore the Fourier Coefficients $\{b_n\}$ are given by:

$$
\begin{equation}
b_k  =  \frac{2}{L}\int\limits_0^Lf(x)\sin\left(\frac{k\pi
x}{L}\right)\, dx.\label{eq:FourierCoefficient}
\end{equation}
$$(ref-fourier-sine-16)

````{prf:example}
:label: example-fourier-sine-0
$$
\begin{eqnarray}
f(x) & = & \left\{
\begin{array}{lll}
2x \quad &0<x<\frac{1}{2}\quad &L=1\\
2(1-x) &\frac{1}{2}<x<1 &\end{array}\right.\nonumber\\
b_n & = & 2\left\{\int\limits_0^{\frac{1}{2}} 2x\sin (n\pi x)\, dx+\int\limits_{\frac{1}{2}}^1 2(1-x)\sin (n\pi x)\, dx\right\}\nonumber\\
& = & 8\frac{\sin (n\pi
/2)}{n^2\pi^2}\qquad\qquad\begin{array}{rcccrcc}
n&= &1&2&3&4&5\\
\sin\left(\frac{n\pi}{2}\right) &&1&0&-1&0&1\end{array}\nonumber\\
\mbox{Therefore } u(x,t) & =&\frac{8}{\pi^2}\sum\limits_{k=0}^\infty
\frac{(-1)^k}{(2k+1)^2}\sin\big[ (2k+1)\pi x\big]
\{\rm\ e\}^{-(2k+1)^2\pi^2t}.
\end{eqnarray}
$$(ref-fourier-sine-17)

- Observe as $t\rightarrow\infty$ $u(x,t)\rightarrow 0$ (all the heat leaks out).
- $\displaystyle u(x,0)=\frac{8}{\pi^2}\sum\limits_{k=0}^\infty\frac{(-1)^k}{(2k+1)^2}\sin\big[ (2k+1)\pi x\big]$.
- $\displaystyle\frac{\pi^2}{8}=\sum\limits_{k=0}^\infty\frac{1}{(2k+1)^2}\qquad$ by letting $\displaystyle x=\frac{1}{2}\Rightarrow f(x)=1$.

```{figure} ../img/fourier/hat_terms.png
:name: hat_terms
:align: center

```
````

````{prf:example}
:label: example-fourier-sine-1
$$
\begin{eqnarray}
f(x) & = & x\qquad 0<x<1\quad L=1\nonumber\\
b_n & = & 2\int\limits_0^1 x\sin (n\pi x)\, dx=-2\frac{\cos (n\pi )}{n\pi }=2\frac{(-1)^{n+1}}{n\pi }\nonumber\\
\mbox{Therefore }u(x,t) & = &
\frac{2}{\pi}\sum\limits_{n=1}^\infty\frac{(-1)^{n+1}}{n}\sin (n\pi
x)\{\rm\ e\}^{-(n\pi )^2t}.
\end{eqnarray}
$$(ref-fourier-sine-18)

- As $t\rightarrow\infty$ $u(x,t)\rightarrow 0$.\\
- $\displaystyle u(x,0)=\frac{2}{\pi }\sum\limits_{n=1}^\infty\frac{(-1)^{n+1}}{n}\sin (n\pi x)$.
- $\begin{array}{ccccl}
  u\left(\displaystyle\frac{1}{2},0\right) &= &\displaystyle\frac{1}{2}&= &\frac{2}{\pi}\sum\limits_{n=1}^\infty \frac{(-1)^{n+1}}{n}\sin (n\pi /2)\\
  &&&= &\displaystyle\frac{2}{\pi}\sum\limits_{k=0}^\infty\frac{(-1)^k}{(2k+1)}\\
  \displaystyle\frac{\pi }{4} &= &&1-\displaystyle\frac{1}{3}+\frac{1}{5} - \ldots
  .\end{array}$

$$
\begin{eqnarray}\begin{array}{ccc}
k&n&\sin\left(\frac{n\pi}{2}\right)\\
0&1&\phantom{-}1\\
 &2&\phantom{-}0\\
1&3&-1\\
 &4&\phantom{-}0\\
2&5&\phantom{-}1\end{array}
\end{eqnarray}
$$(ref-fourier-sine-19)

```{figure} ../img/fourier/x_terms.png
:name: x_terms
:align: center

```
````
