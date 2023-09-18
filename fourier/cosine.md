
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

# Fourier Cosine Series

In this lecture we use separation of variables to solve the heat
equation subject to Neumann boundary conditions. In this case we
reduce the problem to expanding the initial condition function
$f(x)$ in an infinite series of cosine functions - known as the
Fourier Cosine Series.

```{admonition} Key Concepts
Heat Equation; Neumann Boundary
Conditions; separation of variables; Fourier Cosine Series.
```

## The Heat Equation Subject to Homogenous Neumann Boundary Conditions

We consider the heat equation subject to the following initial and
boundary conditions:

```{figure} ../img/fourier/neumann_bar.png
:name: neumann_bar_cosine
:align: center

Consider a conducting bar with thermal conductivity
$\alpha^2$ that has an initial temperature distribution
$u(x,0)=f(x)$ and whose endpoints are insulated
```

$$
\begin{eqnarray}
\mbox{Heat Equation} &:&\ u_{t} =\alpha ^{2}u_{xx},\ \ 0<x<L \label{eq:heat}\\
\mbox{Boundary Conditions}&:& \ \frac{\partial u(0,t)}{\partial x} = 0=\frac{\partial u(L,t)}{\partial x
} \label{eq:BC}\\
\mbox{Initial Condition} &:& \ u(x,0) =f(x) \label{eq:IC}
\end{eqnarray}
$$(ref-fourier-cosine-heat)

### Separation Of Variables - Fourier Sine Series

Consider the heat conduction in an insulated rod whose endpoints are
insulate for all time and within which the initial temperature is
given by $f(x)$ as shown in {numref}`neumann_bar`.

Fourier's Guess:

$$
\begin{eqnarray}
u(x,t) & = & X(x)T(t)\\
u_t & = & X(x)\dot{T}(t)=\alpha^2
u_{xx}=\alpha^2X^{\prime\prime}(x)T(t)\nonumber
\end{eqnarray}
$$(ref-fourier-cosine-1)

$\div\alpha^2 XT$:

$$
\begin{equation}
\frac{X^{\prime\prime}(x)}{X(x)}=\frac{\dot{T}(t)}{\alpha^2
T(t)}=\mbox{ Constant }=-\lambda^2.
\end{equation}
$$(ref-fourier-cosine-2)

Time equation

$$
\begin{eqnarray}\begin{array}{c}
\dot{T}(t)=-\alpha^2\lambda^2 T(t)\quad\displaystyle\frac{dT}{T}=-\alpha^2\lambda^2\, dt\\
\ln|T| =-\alpha^2\lambda^2t+c\\
T(t)=D\{\rm\ e\}^{-\alpha^2\lambda^2t}.\end{array}
\end{eqnarray}
$$(ref-fourier-cosine-3)

**Case I - Spatial equation assuming that $\lambda\ne 0$:**

$$
\begin{eqnarray}\begin{array}{c}
X^{\prime\prime}(x)+\lambda ^2X(x)=0\\
\mbox{Guess}\quad X(x)=\{\rm\ e\}^{rx}\Rightarrow
(r^2+\lambda^2)\{\rm\ e\}^{rx}=0\quad r=\pm\lambda i\end{array}
\end{eqnarray}
$$(ref-fourier-cosine-4)

$$
\begin{eqnarray}\begin{array}{ccl}
X &= &c_1\{\rm\ e\}^{i\lambda x}+c_{2}\{\rm\ e\}^{-i\lambda x}\\
  &= &A\cos\lambda x+B\sin\lambda x \nonumber\\
X'&= &-A\lambda\sin\lambda x+B \lambda\cos\lambda x  \label{eq:Xsolution}\end{array}
\end{eqnarray}
$$(Xsolution)

Now impose the boundary conditions:

$$
\begin{eqnarray}\begin{array}{lcrclcl}
0&=&\frac{\partial u(0,t)}{\partial x}&=&X'(0)T(t)&\Rightarrow& X'(0)=0\\
0&=&\frac{\partial u(L,t)}{\partial x}&=&X'(L)T(t) &\Rightarrow&
X'(L)=0.\end{array}
\end{eqnarray}
$$(ref-fourier-cosine-6)

Now substitute the solution from {eq}`Xsolution` and use
the fact that we have assumed that $\lambda\ne 0$

$$
\begin{eqnarray}\begin{array}{lcrclcl}
0&=&X'(0)&=&-A\lambda . 0+B \lambda&\Rightarrow& B=0\\
0&=&X'(L)&=&-A\lambda sin \lambda L\lambda&\Rightarrow& \lambda_n=
\left(\frac{n\pi}{L}\right)\quad n=1,2,\ldots \end{array}
\end{eqnarray}
$$(ref-fourier-cosine-7)

Therefore for the case $\lambda \ne 0$ we have the countably infinite set of eigenvalues and eigenfunctions

$$
\begin{equation}
\lambda_n =\left(\frac{n\pi}{L}\right)\quad n=1,2,\ldots  \mbox{and }\ X_n(x)= \cos\left(\frac{n\pi x}{L}\right).
\end{equation}
$$(ref-fourier-cosine-8)

**Case II - Spatial equation assuming that $\lambda= 0$:**

In this case the spatial ODE reduces to

$$
\begin{eqnarray}\begin{array}{lcrclcl}
X^{\prime\prime}(x)&=&0\\
\end{array}
\end{eqnarray}
$$(ref-fourier-cosine-9)

which has a general solution

$$
\begin{eqnarray}\begin{array}{lcrclcl}
X(x)&=&A.1+Bx \\
X'(x)&=&B \end{array}
\end{eqnarray}
$$(ref-fourier-cosine-10)

Now imposing the boundary conditions

$$
\begin{eqnarray}\begin{array}{lcrclcl}
0&=&X'(0)&=&B&\Rightarrow& B=0\\
0&=&X'(L)&=&B&\Rightarrow& B=0 \end{array}
\end{eqnarray}
$$(ref-fourier-cosine-11)

The complete set of eigenvalues and eigenfunctions are thus:

$$
\lambda_n =\left(\frac{n\pi}{L}\right)\quad n=0,1,2,\ldots  \mbox{and }\ \ X_0(x) = 1,\ X_n(x)= \cos\left(\frac{n\pi x}{L}\right),\ \quad n=1,2,\ldots.
$$(ref-fourier-cosine-12)

And so,

$$
\quad u_n(x,t)=\{\rm\ e\}^{-\alpha^2{\left(\frac{n\pi}{L}\right)}^2t} \cos\left(\frac{n\pi x}{L}\right)
\quad n=0,1,2,\ldots\nonumber\\
$$(ref-fourier-cosine-13)

are all solutions of $u_t=\alpha^2u_{xx}$ that satisfy the
boundary conditions in {eq}`ref-fourier-cosine-heat`.

Since {eq}`ref-fourier-cosine-heat` is linear, a linear combination of solutions
is again a solution. Thus the most general solution is for the form

$$
\begin{equation}
u(x,t)=A_0+\sum\limits_{n=1}^\infty A_n\cos\left(\frac{n\pi
x}{L}\right)\{\rm\ e\}^{-\alpha^2{\left(\frac{n\pi}{L}\right)}^2t}.
\label{eq:uxtsolution}
\end{equation}
$$(uxtsolution)

What about the initial condition $u(x,0)=f(x)$? If we let $t=0$ in
{eq}`uxtsolution`, then to complete the solution process we are
reduced to determining the coefficients $A_n$ in the series

$$
\begin{equation}
u(x,0)=f(x)=A_0+\sum\limits_{n=1}^\infty A_n\cos\left(\frac{n\pi
x}{L}\right).
\end{equation}
$$(ref-fourier-cosine-15)

As in the last lecture we use the inner product $<.,.>$ to
project $f(x)$ onto the basis functions in the series:

$$
\begin{eqnarray}
& & f(x)=A_0+\sum\limits_{n=1}^\infty A_n\cos\left(\frac{n\pi x}{L}\right)\\
& &\langle f,\cos\left(\frac{k\pi x}{L}\right)\rangle=
\int\limits_0^L f(x)\cos\left(\frac{k\pi x}{L}\right)\,
dx=A_0\int\limits_0^L\cos\left(\frac{k\pi
x}{L}\right)dx+\sum\limits_{n=1}^\infty A_n\int\limits_0^L
\cos\left(\frac{n\pi x}{L}\right)\cos\left(\frac{k\pi x}{L}\right)\,
dx.\label{eq:cosprojection}
\end{eqnarray}
$$(cosprojection)

Recall the identity $\displaystyle\cos (A)\cos B=\frac{1}{2}\left\{\cos
(A-B)+\cos (A+B)\right\}$. Therefore

$$
\begin{eqnarray}
J_{nk} & = & \int\limits_0^L\cos\left(\frac{n\pi x}{L}\right)\cos\left(\frac{k\pi x}{L}\right)\, dx\nonumber\\
& = & \frac{1}{2}\int\limits_0^L\cos (n-k)\frac{\pi x}{L}+\cos (n+k)\frac{\pi x}{L}\, dx\qquad n\not= k\nonumber\\
& = & \frac{1}{2}{\left[\frac{\sin (n-k)\pi x/L}{(n-k)\pi /L}+\frac{\sin (n+k)\pi x/L}{(n+k)\pi /L}\right]}_0^L\nonumber\\
& = & 0\\
J_{nn} & = & \int\limits_0^L\cos^2 \left(\frac{n\pi x}{L}\right)\, dx=\frac{1}{2}\int\limits_0^L 1+\cos\left(\frac{2n\pi x}{L}\right)\, dx\nonumber\\
& = & L/2\nonumber \\
J_{00} & = & \int\limits_0^L 1 dx=L \nonumber
\end{eqnarray}
$$(ref-fourier-cosine-17)

Substituting these integrals into {eq}`cosprojection` we
obtain the following expressions for the Fourier Coefficients $A_k$

$$
\begin{eqnarray}
A_0  &=&  \frac{1}{L}\int\limits_0^L f(x)
dx.\label{eq:FourierCoefficientA0} \\
 A_k  &=&
\frac{2}{L}\int\limits_0^Lf(x)\cos\left(\frac{k\pi x}{L}\right)\,
dx.\label{eq:FourierCoefficientAk}
\end{eqnarray}
$$(FourierCoefficient)

Finally the solution of the initial boundary value problem {eq}`ref-fourier-cosine-heat` is

$$
\begin{equation}
u(x,t)=A_0+\sum\limits_{n=1}^\infty A_n\cos\left(\frac{n\pi
x}{L}\right)\{\rm\ e\}^{-\alpha^2{\left(\frac{n\pi}{L}\right)}^2t}.
\end{equation}
$$(ref-fourier-cosine-19)

where $A_n$ are defined in
{eq}`FourierCoefficient`. We
observe that as $t\rightarrow \infty$  it follows  that $u(x,t)
\rightarrow A_0$, which is just the average value of the initial
heat $f(x)$ distributed in the bar as can be seen from
{eq}`FourierCoefficient`. This is consistent with physical
intuition.

It is sometimes convenient to re-define the Fourier coefficients as follows:

$$
\begin{eqnarray}
a_0  &=&  2 A_0 \nonumber \\
 a_k &=& A_k, \ k=1,2,\dots \nonumber \\
 \mbox{so that the $a_k$ assume the unified form } a_k &=& \frac{2}{L}\int\limits_0^Lf(x)\cos\left(\frac{k\pi x}{L}\right)\,
dx \ \ k=0,1,2,\dots \label{eq:FourierCoefficientsmallak}
\end{eqnarray}
$$(FourierCoefficientsmallak)

In terms of the new coefficients $a_k$ defined in
{eq}`FourierCoefficientsmallak` the Fourier expansion for the
initial condition function $f(x)$ is of the form

$$
\begin{equation}
f(x)=\frac{a_0}{2}+\sum\limits_{n=1}^\infty a_n\cos\left(\frac{n\pi
x}{L}\right) \label{eq:FSCosine}
\end{equation}
$$(FSCosine)

while the solution of the heat equation {eq}`ref-fourier-cosine-heat` is of the form

$$
\begin{equation}
u(x,t)=\frac{a_0}{2}+\sum\limits_{n=1}^\infty
a_n\cos\left(\frac{n\pi
x}{L}\right)\{\rm\ e\}^{-\alpha^2{\left(\frac{n\pi}{L}\right)}^2t}.
\end{equation}
$$(ref-fourier-cosine-22)

````{prf:example}
:label: example-fourier-cosine-0 Fourier Cosine Expansion
Determine the Fourier coefficients $a_k$ for the function

$$
\begin{equation}
f(x) = x, \quad 0<x<1=L
\end{equation}
$$(ref-fourier-cosine-23)

and use the resulting Fourier Cosine expansion to prove the
identity

$$
\frac{\pi ^{2}}{8}=1+\frac{1}{3^{2}}+\frac{1}{5^{2}}+\frac{1}{7^{2}}\ldots +
\frac{1}{\left( 2k+1\right) ^{2}}+\ldots
$$

Solution:

$$
\begin{eqnarray*}
a_{0} &=&2\int_{0}^{1}xdx=2\left[ \frac{x}{2}\right] _{0}^{1}=1 \\
a_{n} &=&2\int_{0}^{1}x\cos (n\pi x)dx=2\frac{(-1)^{n}-1}{n^{2}\pi ^{2}}
=\left\{
\begin{array}{c}
-\frac{4}{n^{2}\pi ^{2}},\ n \ \mbox{odd} \\
0,\ n \ \mbox{even}
\end{array}
\right.
\end{eqnarray*}
$$(ref-fourier-cosine-24)

substituting these expressions for the $a_n$ into {eq}`FSCosine`, we obtain

$$
f(x)=x=\frac{1}{2}-\frac{4}{\pi ^{2}}\sum\limits_{k=0}^\infty
\frac{1}{(2k+1)^{2}}\cos\left({(2k+1)\pi x}\right)
$$(FSCosinex)

To obtain the required identity we set $x=1$ in
and rearrange terms. The partial sums are shown in
{numref}`cosines_close`

```{figure} ../img/fourier/cosines_close.png
:name: cosines_close
:align: center

These figures show the partial sums of the Fourier Cosine Series.
```

In {numref}`cosines_far` we plot the
same graphs but on a larger domain than $[0,L]=[0,1]$.

```{figure} ../img/fourier/cosines_far.png
:name: cosines_far
:align: center

These figures show the partial sums of the Fourier Cosine Series.
```
````
