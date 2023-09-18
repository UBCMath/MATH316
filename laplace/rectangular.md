$$
\newenvironment{remark}{\begin{Remark}\rm}{\end{Remark}}
\newenvironment{exercise}{\begin{Exercise}\rm}{\end{Exercise}}
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

# More Rectangular Domains: Neumann Problems, Mixed BCs, And Semi-Infinite Strip Problems

In this lecture we Proceed with the solution of Laplace's equations
on rectangular domains with Neumann, mixed boundary conditions, and
on regions which comprise a semi-infinite strip.

```{admonition} Key Concepts
Laplace's equation; Rectangular
domains; The Neumann Problem; Mixed BC and semi-infinite strip
problems.
```

Reference Section: Boyce and Di Prima Section 10.8

## The Neumann Problem On A Rectangle - Only Flux Boundary Conditions

````{prf:example}
:label: example-laplace-rectangular-0 The Neumann Problem

```{figure} ../img/laplace/neumann.png
:name: neumann
:align: center

Inhomogeneous Neumann Boundary conditions on a rectangular
domain as prescribed in {eq}`LaplaceNeumannGen`.
```

$$
\begin{eqnarray}
u_{xx}+u_{yy} & = & 0, \qquad 0<x<a\quad 0<y<b\\
u_x(0,y) & = & 0\quad u_x(a,y)=f(y)\\
u_y(x,0) & = & 0=u_y(x,b).
\end{eqnarray}
$$(LaplaceNeumannGen)
````

Let $u(x,y)=X(x)Y(y)$.

$$
\begin{eqnarray}
\frac{X^{\prime\prime}(x)}{X(x)}=
   -\frac{Y^{\prime\prime}(y)}{Y(y)}=\lambda ^2
\end{eqnarray}
$$(ref-laplace-rectangular-1)

$$
\begin{eqnarray}
\left.\begin{array}{l}Y^{\prime\prime}(y)+\lambda^2 Y(y)=0\\
Y^\prime (0)=0=Y^\prime (b)\end{array}\right\}\quad
\begin{array}{lcl}
Y&=&A\cos\lambda y+B\sin\lambda y\\
Y^\prime&=&-A\lambda\sin\lambda y+B\lambda \cos\lambda y\end{array}
\end{eqnarray}
$$(ref-laplace-rectangular-2)

$$
\begin{eqnarray}
Y^\prime (0) & = & \lambda B=0\quad \lambda =0\mbox{ or }B=0.\\
\nonumber\\
Y^\prime (b) & = & -A\lambda\sin\lambda b=0\quad\begin{array}{lcl}
\lambda_n &=&(n\pi /b)\quad n=0,1,\ldots\\
Y_n&=&\cos\left(\displaystyle\frac{n\pi y}{b}\right) ,\quad Y_0=1\end{array}
\end{eqnarray}
$$(ref-laplace-rectangular-3)

$$
\begin{eqnarray}
X_n^{\prime\prime} -\lambda^2 X_n & = & 0\\
X_n^\prime (0) & = & 0
\end{eqnarray}
$$(ref-laplace-rectangular-4)

- $n=0$: $X_0^{\prime\prime}=0$

   $$
   X_0=c_0x+D_0\Rightarrow X_0^\prime
   =c_0\Rightarrow X_0^\prime (0)=c_0=0
   $$

- Choose $D_0=1$: $X_0=1$

   $$
   \begin{eqnarray}
   \begin{array}{rcl}
   n\geq 1\quad X_n&=&c_n\cosh(\lambda_n x)+D_n\sinh(\lambda_n x)\\
   X_n^\prime &=&c_n\lambda \sinh(\lambda_n x)+D_n\lambda\cosh(\lambda_n x)\\
   X_n^\prime (0)&=&\lambda_n D_n=0\end{array}
   \end{eqnarray}
   $$(ref-laplace-rectangular-5)

- Choose $c_n=1$: $X_n=\cosh(\lambda_n x)$

   Thus

   $$
   \begin{eqnarray}
   \left. \begin{array}{lclcl}
   u_n(x,y)&=&X_nY_n&=&\cosh(\lambda_n x)\cos (\lambda_n y)\\
   u_0(x,y)&=&X_0Y_0&=&1\end{array}\right\}\mbox{ satisfy homog.
   BC.\quad}
   \end{eqnarray}
   $$(ref-laplace-rectangular-6)

   Therefore

   $$
   \begin{eqnarray}
   u(x,y)=A_0+\sum\limits_{n=1}^\infty A_n\cosh\left(\frac{n\pi
   x}{b}\right)\cos\left(\frac{n\pi y}{b}\right) .
   \end{eqnarray}
   $$(ref-laplace-rectangular-7)

   Now $f(y)=u_x(a,y)$.

   $$
   \begin{eqnarray}
   u_x(x,y) & = & \sum\limits_{n=1}^\infty
   A_n\left(\frac{n\pi}{b}\right)\sinh
      \left(\frac{n\pi x}{b}\right)\cos\left(\frac{n\pi y}{b}\right)\\
   u_x(a,y) & = & \sum\limits_{n=1}^\infty \left\{
   A_n\left(\frac{n\pi}{b}\right) \sinh
      \left(\frac{n\pi a}{b}\right)\right\} \cos\left(\frac{n\pi y}{b}\right)
      =f(y)\ldots\phantom{\int\int} \label{eqNeumannExpansion}
   \end{eqnarray}
   $$(NeumannExpansion)

   This is like a Fourier Cosine Series for $f(y)$ but without the
   constant term $a_0$.

   Recall

   $$
   \begin{eqnarray}
   f(y)=\frac{a_0}{2}+\sum\limits_{n=1}^\infty a_n\cos\left(\frac{n\pi
   y}{b}\right) ,\; a_n=\frac{2}{b}\int\limits_0^b
   f(y)\cos\left(\frac{n\pi y}{b}\right)\, dy.
   \end{eqnarray}
   $$(ref-laplace-rectangular-9)

   Thus the expansion {eq}`NeumannExpansion` is consistent only if
   $a_0=0$. For this to be true we require that

   $$
   \begin{eqnarray}
   \int\limits_0^b f(y)\, dy=0
   \end{eqnarray}
   $$(ref-laplace-rectangular-10)

if $\displaystyle\int\limits_0^b f(y)\, dy\ne 0$ then there is no solution to
the boundary value problem~1.

```{note}
1. If $\displaystyle\int\limits_0^b f(y)\, dy\ne 0$ there is a {\bf{net flux}}
into the domain through the right hand boundary and, since the other
boundaries are insulated, there can be no steady solution -- the
temperature will continually change with time.

2. If $\displaystyle\int\limits_0^b f(y)\, dy=0$ there is no net flux
through the boundary and a steady state can exist. i.e. It is
possible that $u_{xx}+u_{yy}=u_t=0$. If $\displaystyle\int\limits_0^b f(y)\,
dy=0$ then

   $$
   \begin{eqnarray}
   A_n\left(\frac{n\pi}{b}\right)\sinh\left(\frac{n\pi a}{b}\right)
      = \frac{2}{b}\int\limits_0^b f(y)\cos\left(\frac{n\pi y}{b}\right)\, dy.
   \end{eqnarray}
   $$(ref-laplace-rectangular-11)

   Therefore

   $$
   \begin{eqnarray}
   A_n=\frac{2}{n\pi\sinh\left(\frac{n\pi a}{b}\right)}
   \int\limits_0^b f(y)\cos\left(\frac{n\pi y}{b}\right)\, dy\quad n\geq 1
   \end{eqnarray}
   $$(ref-laplace-rectangular-12)

   and

   $$
   \begin{eqnarray}
   u_\infty (x,y)=A_0+\sum\limits_{n=1}^\infty A_n\cosh
      \left(\frac{n\pi x}{a}\right) \cos\left(\frac{n\pi y}{b}\right)
   \end{eqnarray}
   $$(ref-laplace-rectangular-13)

   where $A_0$ is undetermined. $u(x,y)$ is said to be known up to an
   arbitrary constant.

3. If $u_\infty (x,y)$ is the steady state of a 2D Heat Equation
$u_t=u_{xx}+u_{yy}$ with $u(x,y,0)=u_0(x,y)$ then

   $$
   \begin{eqnarray}
   \int\limits_D u_t\, dx\, dy=\int\limits_D \mathbf{\nabla}\cdot
   \mathbf{\nabla} u\, dx\, dy=\int\limits_{\partial D}\frac{\partial u}{\partial n}\
   ds=0.
   \end{eqnarray}
   $$(ref-laplace-rectangular-14)

   Therefore

   $$
   \begin{eqnarray}
   \frac{\partial}{\partial t}\left(\int\limits_D u\, dx\, dy\right)
   =0\Rightarrow\int\limits_D u\, dx\, dy=\mbox{ const for all time
   }=\int\limits_D u_0(x,y)\, dx\, dy.\phantom{\int}
   \end{eqnarray}
   $$(ref-laplace-rectangular-15)

   Now

   $$
   \begin{eqnarray}
   \int\limits_{ D} u_\infty (x,y)dx dy =A_0 \times \mbox{
   area}(D)=\int\limits_{ D} u_0(x,y)\, dx
   \end{eqnarray}
   $$(ref-laplace-rectangular-16)
   Which is the condition that determines $A_0$.
```

## Rectangular Domains With Mixed BC

````{prf:example}
:label: example-laplace-rectangular-1 Insulating BC along two sides and specified temperatures on the others

```{figure} ../img/laplace/mixed.png
:name: mixed
:align: center

Mixed Boundary conditions on a rectangular domain as
prescribed in {eq}`LaplaceMixed`.
```

$$
\begin{eqnarray}
\Delta u & = &u_{xx}+u_{yy}=0\\
0 & = & u_x(0,y)=u_x(a,y)=u(x,0) \label{eqLaplaceMixed}\\
   &    & u(x,b)=f(x).
\end{eqnarray}
$$(LaplaceMixed)

Let $u(x,y)=X(x)Y(y)$.

$$
\begin{eqnarray}
\frac{X^{\prime\prime}}{X}=-\frac{Y^{\prime\prime}}{Y}=\pm \lambda^2.
\end{eqnarray}
$$(ref-laplace-rectangular-18)

Since we have homogeneous BC on $X^\prime (0)=0=X^\prime (a)$ choose
$-\lambda^2$.

1. $X^{\prime\prime}+\lambda^2 X=0\quad X^{\prime}(0)=0=X^\prime (a)$.

   $$
   \begin{eqnarray}\begin{array}{lcl}
   X(x)&=&A\cos\lambda x+B\sin\lambda x\\
   X^\prime (0)&=&B\lambda =0\Rightarrow B=0\end{array}\begin{array}{lcl}
   X^\prime (x)&=&-A\lambda \sin (\lambda x)+B\lambda \cos (\lambda x)\\
   X^\prime (a)&=&-A\lambda\sin (\lambda a)=0\end{array}
   \end{eqnarray}
   $$(ref-laplace-rectangular-19)

   Therefore

   $$
   \begin{eqnarray}
   \lambda_n =(n\pi /a)\quad n=0,1,2,\ldots\quad
   X_n(x)=\cos\left(\frac{n\pi y}{a}\right)
   \end{eqnarray}
   $$(ref-laplace-rectangular-20)

   are eigenfunctions and eigenvalues.

2. $\lambda_n \ne 0$: $Y^{\prime\prime}-\lambda^2 Y=0$ and $Y(0)=0\Rightarrow Y_n(y)=A\sinh\displaystyle\left(\frac{n\pi y}{a}\right)$ $n\ne 0$. Thus

   $$
   \begin{eqnarray}
   u_n(x,y)=\cos\left(\frac{n\pi x}{a}\right) \sinh\left(\frac{n\pi
   y}{a}\right)
   \end{eqnarray}
   $$(ref-laplace-rectangular-21)

   satisfy homogeneous BC.

${\la_0 =0}$: In this case the ODE for $Y_0$ is:

$$
\begin{eqnarray}
Y_0^{\prime\prime}&=&0\Rightarrow Y(y)=c_1y+c_2\\
Y_0(0)&=&c_2=0\Rightarrow Y_0(y)=y
\end{eqnarray}
$$(ref-laplace-rectangular-22)

and $u_0(x,y)=y\cdot 1$ satisfies the homogeneous BC.

Therefore

$$
\begin{eqnarray}
u(x,y) & = & c_0y+\sum\limits_{n=1}^\infty c_n\sinh\left(\frac{n\pi y}{a}\right)\cos\left(\frac{n\pi x}{a}\right)\\
u(x,b) & = & \frac{(2c_0b)}{2}+\sum\limits_{n=1}^\infty c_n\sinh\left(\frac{n\pi b}{a}\right)\cos\left(\frac{n\pi x}{a}\right) =f(x)\\
(2c_0b) & = & \frac{2}{a}\int\limits_0^a f(x)\, dx;\
\begin{array}{l}
c_n\sinh\left(\frac{n\pi b}{a}\right) =\frac{2}{a}\int\limits_0^a f(x)\cos \left(\frac{n\pi x}{a}\right)\, dx\phantom{\int}\end{array}\\
c_0 & = & \frac{1}{ab}\int\limits_0^a f(x)\, dx;\ \begin{array}{l}
c_n=\frac{2}{a\sinh\left(\frac{n\pi b}{a}\right)}\int\limits_0^a
f(x)\cos\left(\frac{n\pi x}{a}\right)\, dx\end{array}\\
u(x,y) & = & c_0y+\sum\limits_{n=1}^\infty c_n\sinh\left(\frac{n\pi
y}{a}\right)\cos\left(\frac{n\pi x}{a}\right) .
\end{eqnarray}
$$(ref-laplace-rectangular-23)
````

## Semi-Infinite Strip Problems

````{prf:example}
:label: example-laplace-rectangular-2 A Semi-Infinite Strip with Specified Temperatures

```{figure} ../img/laplace/strip_homo.png
:name: strip
:align: center

Diriclet Boundary conditions on a semi-infinite strip as
prescribed in (\ref{semiHomog}).
```

$$
\begin{eqnarray}
u_{xx}+u_{yy} & = & 0\qquad 0<x<a,\quad 0<y<\infty\\
u(0,y) & = & 0=u(a,y)\\
u(x,0) & = & f(x)\qquad u(x,y)\rightarrow 0\mbox{ as }y\rightarrow\infty
\label{semiHomog}
\end{eqnarray}
$$(ref-laplace-rectangular-24)

Let $u(x,t)=X(x)T(t)$ and plug into (1a?):

$$
\begin{eqnarray}
\frac{X^{\prime\prime}(x)}{X(x)}= -
\frac{Y^{\prime\prime}(y)}{Y(y)}=-\lambda^2 \mbox{ since we have
homogeneous BC on $X$}.
\end{eqnarray}
$$(ref-laplace-rectangular-25)

1.   
   $$
   \begin{eqnarray}\left.\begin{array}{l}
   X^{\prime\prime}+\lambda^2 X=0\\
   X(0)=0=X(a)\end{array}\right\}\quad\begin{array}{lcl}
   \lambda_n &=&n\pi /a\quad n=1,2, \ldots\\
   X_n &=&\sin\left(\displaystyle\frac{n\pi x}{a}\right)\end{array}
   \end{eqnarray}
   $$(ref-laplace-rectangular-26)

2. $Y^{\prime\prime} -\lambda^2 Y=0\quad Y(y)=A\{\rm\ e\}^{-\lambda y}+B\{\rm\ e\}^{\lambda y}$. Since $u(x,y)\rightarrow 0$ as $y\rightarrow\infty$ we require $B=0$. Therefore

   $$
   \begin{eqnarray}
   u_n(x,y)=\{\rm\ e\}^{-\lambda_n y}\sin\left(\frac{n\pi x}{a}\right)
   \end{eqnarray}
   $$(ref-laplace-rectangular-27)
   satisfy the homogeneous BC and the BC at $\infty$. Thus

   $$
   \begin{eqnarray}
   u(x,y)=\sum\limits_{n=1}^\infty
   c_n\{\rm\ e\}^{-\left(\frac{n\pi}{a}\right) y}\sin\left(\frac{n\pi
   x}{a}\right) .
   \end{eqnarray}
   $$(ref-laplace-rectangular-28)

   $$
   \begin{eqnarray}
   f(x)=u(x,0)=\sum\limits_{n=1}^\infty c_n\sin\left(\frac{n\pi
   x}{a}\right)\Rightarrow c_n=\frac{2}{a}\int\limits_0^a
   f(x)\sin\left(\frac{n\pi x}{a}\right)\, dx.\phantom{\int}
   \end{eqnarray}
   $$(ref-laplace-rectangular-29)
````

````{prf:example}
:label: example-laplace-rectangular-3 Semi-Infinite Strip with Inhomogeneous BC

```{figure} ../img/laplace/strip_inhomo.png
:name: strip_inhomo

Diriclet Boundary conditions on a semi-infinite strip as
prescribed in (\ref{semiInHomog}).
```

$$
\begin{eqnarray}
u_{xx}+u_{yy} & = & 0\qquad 0<x<a,\quad 0<y<\infty\\
u(0,y) & = & A,\quad B=u(a,y)\\
u(x,0) & = & f(x)\qquad u(x,y)\rightarrow 0\mbox{ as }y\rightarrow\infty
\label{semiInHomog}
\end{eqnarray}
$$(ref-laplace-rectangular-30)

Look for a function $v(x)$ for which $v^{\prime\prime}=0$ and which
satisfies the inhomogeneous BC.

$v=\alpha x+\beta\quad v(0)=A=\beta\quad v(a)=\alpha a+A=B$

Therefore $v(x)=\displaystyle\left(\frac{B-A}{a}\right) x+A$.

Now let $u(x,y)=v(x)+w(x,y)$.

$$
\begin{eqnarray}
0 & = & u_{xx}+u_{yy}=v_{xx}\!\!\!\!\!\!\!\nearrow +w_{xx}+v_{yy}\!\!\!\!\!\!\!\nearrow +w_{yy}\Rightarrow\Delta w=0\\
A & = & u(0,y)=v(0)+w(0,y)\Rightarrow w(0,y)=0\\
B & = & u(a,y)=v(a)+w(a,y)\Rightarrow w(a,y)=0\\
f(x) & = & u(x,0)=v(x)+w(x,0)\Rightarrow w(x,0)=f(x)-v(x).
\end{eqnarray}
$$(ref-laplace-rectangular-31)

Thus $w$ satisfies the same BVP as does $u$ in Eg.~3 above.

Therefore

$$
\begin{eqnarray}
u(x,y)=(B-A)(x/a)+A+\sum\limits_{n=1}^\infty
d_n\{\rm\ e\}^{-\left(\frac{n\pi}{a}\right) y}\sin\left(\frac{n\pi
x}{a}\right)
\end{eqnarray}
$$(ref-laplace-rectangular-32)

where

$$
\begin{eqnarray}
d_n=\frac{2}{a}\int\limits_0^a\left\{
f(x)-v(x)\right\}\sin\left(\frac{n\pi x}{a}\right)\, dx.
\end{eqnarray}
$$(ref-laplace-rectangular-33)
````
