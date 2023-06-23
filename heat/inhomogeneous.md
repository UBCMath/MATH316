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

# Heat Conduction Problems With Time-Independent Inhomogeneous Boundary Conditions

In this lecture we consider heat conduction problems with
inhomogeneous boundary conditions. To determine a solution we
exploit the linearity of the problem, which guarantees that linear
combinations of solutions are again solutions. In particular,
(excuse the pun) we first determine a well chosen particular
solution, known as the steady-state solution that can be used
to remove the inhomogeneous boundary conditions. This reduces the
problem to one of solving the same boundary value problem but with
homogeneous boundary conditions and an augmented initial
condition. Although the steady state solution is a natural choice in
this case, the choice of particular solution, as always, is by
no means unique. To solve the homogeneous boundary value problems we
demonstrate two distinct methods: Method I: comprises the more
elementary method of separation of variables; while Method II
introduces the more generally applicable method of
eigenfunction expansions.

```{admonition} Key Concepts
Inhomogeneous Boundary Conditions,
Particular Solutions, Steady state Solutions; Separation of
variables, Eigenvalues and Eigenfunctions, Method of Eigenfunction
Expansions.
```

Reference Sections: Boyce and Di Prima Sections 10.5, 10.6, 11.2, and 11.3.

## Heat Conduction Problems with Inhomogeneous Boundary Conditions

### A Summary Of Eigenvalue Boundary Value Problems And Their Eigenvalues And Eigenfunctions

Thus far we have discussed five fundamental Eigenvalue problems:

1. The Dirichlet Problem
2. The Neumann Problem
3. Periodic Boundary Conditions
4. Two types of Mixed Boundary Value Problems

Since these Eigenvalue problems will recur throughout the remainder of the
course, for convenient reference we list the boundary value ODEs and
the corresponding eigenvalues and eigenfunctions. We also plot the
first few eigenfunctions, which can be seen to satisfy the
prescribed boundary conditions. You will see that it is sometimes
convenient to remember these so that we do not need to resort to the
method of separation of variables in order to derive the appropriate
eigenvalues and eigenfunctions for a particular problem. Recognizing
the appropriate eigenvalues and eigenfunctions for a given problem
gives rise to the so-called method of eigenfunction
expansions, which we introduce in it simplest form in this lecture.
There is, however, a word of warning. If the boundary value problem
is new, in that it includes extra terms that cannot be grouped with
the time variables in the separation process, it is necessary to
consider the new eigenvalue problem in order to determine the
appropriate eigenfunctions and eigenvalues.

#### The Dirichlet Problem

```{figure} ../img/heat/dirichlet_modes.png
:name: dirichlet_modes
:alt: dirichlet_modes
:align: center

Dirichlet BCs, e.g. with ice on both sides.
```

$$
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X(0)=0=X(L)
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\lambda _{n}=\frac{n\pi }{L},\ n=1,2,\ldots  \\
X_{n}(x)=\sin \left( \frac{n\pi x}{L}\right)
\end{array}
\right.
$$(ref-heat-inhomogeneous-0)

#### The Neumann Problem

```{figure} ../img/heat/neumann_modes.png
:name: neumann_modes
:alt: neumann_modes
:align: center

Neumann BCs, e.g. with insulation on both sides.
```

$$
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X^{\prime }(0)=0=X^{\prime }(L)
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\lambda _{n}=\frac{n\pi }{L},\ n=0,1,2,\ldots  \\
X_{n}(x)=\cos \left( \frac{n\pi x}{L}\right)
\end{array}
\right.
$$(ref-heat-inhomogeneous-1)

#### The Periodic Boundary Value Problem

```{figure} ../img/heat/periodic_modes.png
:name: periodic_modes
:alt: periodic_modes
:align: center

Periodic BCs, e.g. a closed ring.
```

$$
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X(-L)=X(L) \\
X^{\prime }(-L)=X^{\prime }(L)
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\lambda _{n}=\frac{n\pi }{L},\ n=0,1,2,\ldots  \\
X_{n}(x)\in \left\{ 1,\cos \left( \frac{n\pi x}{L}\right) ,\sin \left( \frac{
n\pi x}{L}\right) \right\}
\end{array}
\right.
$$(ref-heat-inhomogeneous-2)

#### Mixed Boundary Value Problem - Type A

```{figure} ../img/heat/mixed_modes_A.png
:name: mixed_modes_A
:alt: mixed_modes_A
:align: center

Mixed BCs Type A, e.g. ice left and insulation right.
```

$$
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X(0)=0=X^{\prime }(L)
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\lambda _{k}=\frac{(2k+1)\pi }{2L},\ k=0,1,2,\ldots  \\
X_{n}(x)=\sin \left( \frac{(2k+1)\pi }{2L}x\right)
\end{array}
\right.
$$(ref-heat-inhomogeneous-3)

#### Mixed Boundary Value Problem - Type B

```{figure} ../img/heat/mixed_modes_B.png
:name: mixed_modes_B
:alt: mixed_modes_B
:align: center

Mixed BCs Type B, e.g. insulation left and ice right.
```

$$
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X^{\prime }(0)=0=X(L)
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\lambda _{k}=\frac{(2k+1)\pi }{2L},\ k=0,1,2,\ldots  \\
X_{n}(x)=\cos \left( \frac{(2k+1)\pi }{2L}x\right)
\end{array}
\right.
$$(ref-heat-inhomogeneous-4)

### Specified Constant Temperatures

````{prf:example}
:label: example-heat-inhomogeneous-0
$$
u_t = \alpha^2 u_{xx}\qquad 0<x<L,\quad t>0 \label{eg16eq1}
$$(eg16eq1)
$$
\mbox{BC: } u(0,t) = u_0\qquad u(L,t)=u_1\quad u_0,u_1\mbox{ constants} \label{eg16eq2}
$$(eg16eq2)

$$
\mbox{IC: } u(x,0) = g(x). \label{eg16eq3}
$$(eg16eq3)

```{figure} ../img/heat/inhomogeneous_dirichlet.png
:name: inhomogeneous_dirichlet
:alt: inhomogeneous_dirichlet
:align: center

Initial, transient, and steady solutions to the heat
conduction problem {eq}`eg16eq3` with inhomogeneous Dirichlet BC
```

__Subtracting out a particular solution:__

Firstly consider the steady-state solution (i.e., when $u_t=0$)
which we denote by $u_{\infty} (x)$. In this case {eq}`eg16eq1`
becomes

$$
\begin{eqnarray*}
\alpha^2 u_\infty ^{\prime\prime}(x) & = & 0\Rightarrow u_\infty (x)=A_0x+B_0\\
u_\infty (0) & = & B_0=u_0\quad u_\infty
(L)=A_0L+u_0=u_1\Rightarrow\quad
\begin{array}{c}
u_\infty (x)=\left(\frac{u_1-u_0}{L}\right) x+u_0\\
\mbox{steady state solution}\end{array}.
\end{eqnarray*}
$$(ref-heat-inhomogeneous-6)

Let $u(x,t)=u_\infty (x)+v(x,t)$. Substitute into {eq}`eg16eq1`

$$
\begin{equation}
u_t= {\big( u_\infty (x)+v(x,t)\big)}_t =\alpha^2 {\big( u_\infty
(x)+v(x,t)\big)}_{xx} \Rightarrow v_t =\alpha^2 v_{xx}
\end{equation}
$$(ref-heat-inhomogeneous-7)

since ${\big( u_\infty (x)\big)}_{xx}=0$. Substitute into {eq}`eg16eq2`

$$
\begin{eqnarray*}\begin{array}{lclclclclcl}
u(0,t) &= &u_0 &= &u_\infty (0)+v(0,t) &= &u_0+v(0,t) &\Rightarrow &v(0,t) &= &0\\
u(L,t) &= &u_1 &= &u_\infty (L)+v(L,t) &= &u_1+v(L,t) &\Rightarrow
&v(L,t) &= &0.\end{array}
\end{eqnarray*}
$$(ref-heat-inhomogeneous-8)

Substitute into {eq}`eg16eq3`

$$
\begin{eqnarray*}
u(x,0)=g(x)=u_\infty (x)+v(x,0)\Rightarrow v(x,0)=g(x)-u_\infty (x).
\end{eqnarray*}
$$(ref-heat-inhomogeneous-9)

Thus we have to solve a new problem for $v$ which has homogenous BC:

$$
\begin{eqnarray}
\begin{array}{lcl}
v_t & = &\alpha ^2 v_{xx}\\
v(0,t)& = & 0=v(L,t)\\
v(x,0)& = &g(x)-u_\infty (x).\end{array}
\end{eqnarray}
$$(ref-heat-inhomogeneous-10)

__Method I: Solving the Homogeneous Problem using Separation of Variables__

Separate variables:
$v(x,t)=X (x)T(t)$.

$$
\begin{eqnarray*}
\frac{\dot{T}(t)}{\alpha^2 T(t)} & = & \frac{X^{\prime\prime}(x)}{X(x)}= -\lambda^2 =\mbox{ const}\\
T(t) & = & c\{\rm\ e\}^{-\lambda^2\alpha^2 t}
\end{eqnarray*}
$$(ref-heat-inhomogeneous-11)

$$
\begin{eqnarray*}
X^{\prime\prime}+\lambda^2 X=0\quad X(0)=0=X(L)\Rightarrow
X_n(x)=\sin\left(\frac{n\pi x}{L}\right),\  \lambda_n
=\frac{n\pi}{L},\ n=1,\ldots
\end{eqnarray*}
$$(ref-heat-inhomogeneous-12)

$$
\begin{eqnarray*}
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X(0)=0=X(L)
\end{array}
\right\}  &\Longrightarrow &\left\{
\begin{array}{c}
\lambda _{n}=\frac{n\pi }{L},n=1,2,\ldots  \\
X_{n}(x)=\sin \left( \frac{n\pi x}{L}\right)
\end{array}
\right.
\end{eqnarray*}
$$(ref-heat-inhomogeneous-13)

Thus

$$
\begin{eqnarray}
v(x,t) & = & \sum\limits_{n=1}^\infty b_n\{\rm\ e\}^{-\alpha^2{\left(\frac{n\pi}{L}\right)}^2 t} \sin\left(\frac{n\pi x}{L}\right)\\
v(x,0) & = & g(x)-u_\infty (x)=\sum\limits_{n=1}^\infty
b_n\sin\left(\frac{n\pi x}{L}\right)\Rightarrow
b_n=\frac{2}{L}\int\limits_0^L\left\{ g(x)-u_\infty
(x)\right\}\sin\left(\frac{n\pi x}{L}\right)\, dx\nonumber
\end{eqnarray}
$$(ref-heat-inhomogeneous-14)

Thus the solution to the inhomogeneous problem is:

$$
\begin{eqnarray}
u(x,t) & = & u_\infty (x)+v(x,t)\\
& = & u_0+\left(\frac{u_1-u_0}{L}\right) x+\sum\limits_{n=1}^\infty
b_n
   \{\rm\ e\}^{-\alpha^2 {\left(\frac{n\pi}{L}\right)}^2 t}\sin\left(\frac{n\pi x}{L}\right) \label{eg161eqres}\\
\mbox{where}& & \nonumber \\
b_n & = & \frac{2}{L}\int\limits_0^L\left\{ g(x)-u_\infty
(x)\right\}\sin\left(\frac{n\pi x}{L}\right)\, dx.
\end{eqnarray}
$$(eg161eqres)

__Method II: Solving the Homogeneous Problem using Eigenfunction Expansions__

In order to solve the boundary value problem
{eq}`eg16eq1`-{eq}`eg16eq3` we could recognize that
$\displaystyle{\left\{\sin\left(\frac{n\pi
x}{L}\right)\right\}}_{n=1}^\infty$ are eigenfunctions of the
spatial operator:

$$
\begin{eqnarray}
-\frac{\partial^2}{\partial x^2}
\end{eqnarray}
$$(ref-heat-inhomogeneous-16)

along with the homogeneous Dirichlet BC $v(0,t)=0=v(L,t)$. We
therefore assume an eigenfunction expansion of the form:

$$
\begin{eqnarray}
v(x,t) & = &\sum\limits_{n=1}^\infty \hat{v}_n(t)\sin\left(\frac{n\pi x}{L}\right)\\
\frac{\partial v}{\partial t} & = &\sum\limits_{n=1}^\infty
\dot{\hat{v}}_n(t)\sin\left(\frac{n\pi x}{L}\right)
  \quad\mbox{and}\quad \\
  \frac{\partial^2 v}{\partial x^2}&=&-\sum\limits_{n-1}^\infty \hat{v}_n(t)
  {\left(\frac{n\pi}{L}\right)}^2 \sin\left(\frac{n\pi x}{L}\right)\\
v_t & = & \alpha^2 v_{xx}\Rightarrow\sum\limits_{n=1}^\infty\left\{
  \dot{\hat{v}}_n(t)+\alpha^2
  {\left(\frac{n\pi}{L}\right)}^2{\hat{v}}_n(t)\right\}
  \sin\left(\frac{n\pi x}{L}\right) =0.\nonumber\\
& &
\end{eqnarray}
$$(ref-heat-inhomogeneous-17)

Therefore

$$
\begin{eqnarray}
\dot{\hat{v}}_n(t) & = & -\alpha^2
{\left(\frac{n\pi}{L}\right)}^2\hat{v}_n(t)\quad
   \mbox{A simple ODE for $\hat{v}_n(t)$:}\\
\Rightarrow \hat{v}_n(t) & = & \hat{v}_n(0)\{\rm\ e\}^{-\alpha^2
{\left(\frac{n\pi}{L}\right)}^2 t}.
\end{eqnarray}
$$(ref-heat-inhomogeneous-18)

Therefore

$$
\begin{eqnarray}
v(x,t) & = &\sum\limits_{n=1}^\infty\hat{v}_n(0)\{\rm\ e\}^{-\alpha^2 {\left(\frac{n\pi}{L}\right)}^2 t}\sin\left(\frac{n\pi x}{L}\right)\\
v(x,0) & = &\sum\limits_{n=1}^\infty \hat{v}_n(0)\sin\left(\frac{n\pi x}{L}\right) = g(x)-u_\infty (x)\\
\hat{v}_n(0) & = & \frac{2}{L}\int\limits_0^{L}\left\{
  g(x)-u_\infty (x)\right\} \sin\left(\frac{n\pi x}{L}\right)\, dx
\end{eqnarray}
$$(ref-heat-inhomogeneous-19)

which is the same solution as that in {eq}`eg161eqres` above.
````

````{prf:example}
:label: example-heat-inhomogeneous-1
$$
\begin{eqnarray}
u_t & = & \alpha^2 u_{xx}\qquad 0<x<L,\quad t>0\\
\mbox{BC: }u(0,t) & = & u_0\quad u_x(L,t)=0\\
\mbox{IC: }u(x,0) & = & g(x).
\end{eqnarray}
$$(ref-heat-inhomogeneous-20)

```{figure} ../img/heat/inhomogeneous_mixed.png
:name: inhomogeneous_mixed
:alt: inhomogeneous_mixed
:align: center

Initial, transient, and steady solutions to the heat
conduction problem {eq}`eg16eq3` with inhomogeneous Mixed BC
```

Look for a steady solution: $u_\infty^{\prime\prime}(x)=0$.

$$
\begin{eqnarray}
u_\infty (x)=Ax+B\quad u_\infty (0)=B=u_0\quad u_\infty ^\prime
(x)=A=0
\end{eqnarray}
$$(ref-heat-inhomogeneous-21)

Therefore

$$
\begin{eqnarray}
u_\infty (x)=u_0.
\end{eqnarray}
$$(ref-heat-inhomogeneous-22)

Let $u(x,t)=u_\infty (x)+v(x,t)=u_0+v(x,t)$.

$$
\begin{eqnarray}
u_t&=&\alpha^2 u_{xx} \Rightarrow   v_t=\alpha^2 v_{xx}\\
u(0,t)&=& u_0 \Rightarrow   u_0 = u_0+v(0,t)\Rightarrow v(0,t)=0\\
u_x(L,t)&=& 0 \Rightarrow   0=v_x(L,t)\Rightarrow v_x(L,t)=0 \\
u(x,0)&=&u_0+v(x,0)=g(x) \Rightarrow   v(x,0)=g(x)-u_0
\end{eqnarray}
$$(ref-heat-inhomogeneous-23)

Thus $v(x,t)$ satisfies

$$
\begin{eqnarray}
v_t&=&\alpha^2 v_{xx} \label{eg162veq}\\
v(0,t)&= &0=v_x(L,t)\\
v(x,0)&=&g(x)-u_0.
\end{eqnarray}
$$(eg162veq)

Further,

$$
u(x,0)=u_0+v(x,0)=g(x)\Rightarrow v(x,0)=g(x)-u_0.
$$

We now need a solution $v(x,t)=X(x)T(t)$ to {eq}`eg162veq`:

$$
\begin{eqnarray}
\frac{\dot{T}(t)}{\alpha^2 T(t)} & = &\frac{X^{\prime\prime}(x)}{X(x)}=-\lambda^2\nonumber\\
\dot{T}(t) & = & -\lambda^2 \alpha^2 T(t)\Rightarrow T(t)=c\{\rm\ e\}^{-\lambda^2 \alpha^2 t}\\
X^{\prime\prime}+\lambda^2 X & = & 0;\quad X(0)=0=X^\prime (L)\\
\Rightarrow X(x) & = & A\cos (\lambda x)+B\sin (\lambda x)\\
X^\prime (x)
  & = & -A\lambda\sin (\lambda x)+B\lambda\cos (dx)\nonumber\\
X(0) & = & A=0\\
\mbox{Therefore }X^\prime (L) & = & B\lambda\cos (\lambda
L)=0\Rightarrow
  \lambda_k = (2k-1)\frac{\pi}{2L},\ k=1,2,3,\ldots\\
\mbox{or } \lambda & = & 0\quad\mbox{which yields the trivial solution.}
\end{eqnarray}
$$(ref-heat-inhomogeneous-25)

Therefore

$$
\begin{eqnarray}
v(x,t) & = & \sum\limits_{k=1}^\infty b_k\{\rm\ e\}^{-\lambda_k^2 \alpha^2 t}\sin\left(\frac{(2k-1)}{2L}\pi x\right)\\
v(x,0) & = & \sum\limits_{k=1}^\infty b_k\sin\left(\frac{(2k-1)}{2L}\pi x\right) =g(x)-u_0\\
\Rightarrow b_n & = & \frac{2}{L}\int\limits_0^L\left\{
g(x)-u_0\right\} \sin\left(\frac{(2k-1)}{2L}\pi x\right)\, dx.
\end{eqnarray}
$$(ref-heat-inhomogeneous-26)

Returning to $u(x,t)=u_0+v(x,t)$:

$$
\begin{eqnarray}
u(x,t)=u_0+\sum\limits_{k=1}^\infty b_k\{\rm\ e\}^{-\alpha^2 \lambda_k^2
t}\sin\left(\frac{(2k-1)}{2L}\pi x\right) .
\end{eqnarray}
$$(ref-heat-inhomogeneous-27)
````

### Heat Conduction With Some Heat Loss And Inhomogeneous Boundary Conditions

````{prf:example} Heat Equation with some Heat Loss
:label: example-heat-inhomogeneous-2

$$
\begin{eqnarray}
u_t & = & \alpha ^2 u_{xx}-u\qquad 0<x<L,\quad t>0\\
\mbox{BC: } u(0,t) & = & 0\quad u(L,t)=u_1\\
\mbox{IC: }u(x,0) & = & g(x).
\end{eqnarray}
$$(ref-heat-inhomogeneous-28)

```{figure} ../img/heat/heat_loss.png
:name: heat_loss
:alt: heat_loss
:align: center

Bar subject to heat loss all along its length with
inhomogeneous Mixed BC
```

Look for the steady state solution $u_\infty (x)$:

$$
\begin{eqnarray}\begin{array}{l}
\alpha^2 u_\infty ^{\prime\prime}-u_\infty =0\\
\displaystyle u_\infty (x)=A\cosh\left(\frac{x}{\alpha}\right) +B\sinh\left(\frac{x}{\alpha }\right)\\
\displaystyle u_\infty (0)=A=0\quad u_\infty (L)=B\sin
h\left(\frac{L}{\alpha}\right) =u_1\quad B=\displaystyle\frac{u_1}{\sin
h\left(\frac{L}{\alpha}\right)} .\end{array}
\end{eqnarray}
$$(ref-heat-inhomogeneous-29)

Therefore

$$
\begin{eqnarray}
u_\infty (x)=u_1\sinh\left.\left(\frac{x}{\alpha}\right)\right/
   \sinh\left(\frac{L}{\alpha }\right) .
\end{eqnarray}
$$(ref-heat-inhomogeneous-30)

```{figure} ../img/heat/steady_state.png
:name: steady_state
:alt: steady_state
:align: center

The steady-state solution.
```

Now let $u(x,t)=u_\infty (x)+v(x,t)$.

$$
\begin{eqnarray}
& &\begin{array}{lc}
u_t=\alpha^2 u_{xx}-u&\Rightarrow\\
u(0,t)=0&\Rightarrow\\
u(L,t)=u_1&\Rightarrow\\
u(x,0)=g(x)&\Rightarrow\end{array}\quad\begin{array}{l}
v_t=\alpha^2 v_{xx}-v\\
0_{\phantom{t}}=u_\infty (0)+v(0,t)\\
u_1=u_\infty (L)+v(L,t)=u_1+v(L,t)\\
u_\infty (x)+v(x,0)=g(x)\end{array}\quad\begin{array}{c}
\\
 \nonumber\\
\\
\\
\end{array}\\
& \Rightarrow &\qquad\left\{\begin{array}{rcl}
v_t&=&\alpha^2 v_{xx}-v\\
v(0,t)&=&0\\
v(L,t)&=&0\\
v(x,0)&=&g(x)-u_\infty (x).\end{array}\right. \label{eg154veq}
\end{eqnarray}
$$(eg154veq)

To solve {eq}`eg154veq` we separate variables $v(x,t)=X(x)T(t)$.
Therefore

$$
\begin{eqnarray}
\frac{\dot{T}(t)}{T(t)}=\frac{\alpha^2
X^{\prime\prime}}{X}-1\Rightarrow
   \frac{1}{\alpha^2}\left(\frac{\dot{T}(t)}{T(t)}+1\right)
   = \frac{X^{\prime\prime}(x)}{X(x)}=-\lambda ^2.
\end{eqnarray}
$$(ref-heat-inhomogeneous-32)

Therefore

$$
\begin{eqnarray}
\dot{T}(t) & = & -(\lambda^2 \alpha^2+1)T(t)\Rightarrow T(t)=c\{\rm\ e\}^{-(1+\lambda^2 \alpha^2 )t}\\
X^{\prime\prime} +\lambda^2 X & = & 0\Rightarrow X(x)=A\cos\lambda x+B\sin\lambda x\\
\Rightarrow X(0) & = & 0\Rightarrow A=0\quad X(L)=B\sin (\lambda
L)=0\Rightarrow\lambda_n =\left(\frac{n\pi}{L}\right)\nonumber\\
& &\hspace{2.2in}n=1,2,\ldots .
\end{eqnarray}
$$(ref-heat-inhomogeneous-33)

Therefore

$$
\begin{eqnarray}
v(x,t) & = & \sum\limits_{n=1}^\infty b_n\{\rm\ e\}^{-(1+\lambda_n^2 \alpha^2) t}\sin\left(\frac{\pi x}{L}\right)\\
v(x,0) & = & g(x)-u_\infty (x)=\sum\limits_{n=1}^\infty
b_n\sin\left(\frac{n\pi x}{L}\right)\Rightarrow
b_n\\
& &\qquad =\frac{2}{L}\int\limits_0^L\left\{ g(x)-u_\infty
(x)\right\}\sin\left(\frac{n\pi x}{L}\right)\, dx.
\end{eqnarray}
$$(ref-heat-inhomogeneous-34)

Therefore

$$
\begin{eqnarray}
u(x,t)=u_1\sinh\left. \left(\frac{x}{\alpha}\right) \right/ \sin
h\left(\frac{L}{\alpha}\right) +\sum\limits_{n=1}^\infty b_n
   \{\rm\ e\}^{-(1+\lambda_n^2 \alpha^2) t}\sin\left(\frac{n\pi x}{L}\right) .
\end{eqnarray}
$$(ref-heat-inhomogeneous-35)

```{note}
The $-u$ term in the PDE is responsible for the $\{\rm\ e\}^{-t}$ factor in the
solution.
```
````

### Heat Conduction with Inhomogeneous Neumann Boundary Conditions

````{prf:example} Inhomogeneous Neumann BC
:label: example-heat-inhomogeneous-3

$$
u_t = \alpha^2 u_{xx}\qquad 0<x<L,\quad t>0\label{NeumannIhomog1}\\
$$(NeumannIhomog1)

$$
\mbox{BC: }u_x(0,t) = A\quad u_x(L,t)=B
$$(NeumannIhomog2)

$$
\mbox{IC: }u(x,0) = g(x) \label{NeumannIhomog3}.
$$(NeumannIhomog3)

```{figure} ../img/heat/inhomogeneous_neumann.png
:name: inhomogeneous_neumann
:alt: inhomogeneous_neumann
:align: center

Initial, transient, and steady solutions to the heat
conduction problem {eq}`NeumannIhomog1`-{eq}`NeumannIhomog3`
with inhomogeneous Neumann BC
```

Try for a steady solution: $u_\infty^{\prime\prime}(x)=0$, $u_\infty (x)=\alpha
x+\beta$, $u_x=\alpha$ but then we cannot match both BC unless
$A=B=\alpha$. This means that if we are pumping and removing heat
from the rod at different rates then the temperature does not reach
a steady state.

Instead of subtracting off a steady solution we subtract a particular solution which depends on $x$ and $t$ of the form:

$$
\begin{eqnarray}
w(x,t) & = & ax^2+bx+ct\\
w_t & = & c=\alpha^2 w_{xx}=2\alpha ^2a\Rightarrow c=2\alpha^2 a.
\end{eqnarray}
$$(ref-heat-inhomogeneous-37)

Then

$$
\begin{eqnarray}
w(x,t)=ax^2+bx+2\alpha^2 at
\end{eqnarray}
$$(ref-heat-inhomogeneous-38)

solves the heat equation.

Now we determine the constants $a$ and $b$ so that $w(x,t)$
satisfies the inhomogeneous BC:

$$
\begin{eqnarray}
w_x=2ax+b:\quad w_x(0,t)=b=A,\quad w_x(L,t)=2aL+A=B.
\end{eqnarray}
$$(ref-heat-inhomogeneous-39)

Therefore $a=(B-A)/2L$. Therefore

$$
\begin{eqnarray}
w(x,t)=\frac{(B-A)}{2L}x^2+Ax+\alpha ^2\left(\frac{B-A}{L}\right) t.
\end{eqnarray}
$$(ref-heat-inhomogeneous-40)

Now let $u(x,t)=w(x,t)+v(x,t)$

$$
\begin{eqnarray}
\begin{array}{lclcl}
u_t&=&w_t+v_t=\alpha^2 (w_{xx}+v_{xx}) &\Rightarrow &v_t=\alpha^2 v_{xx}\\
u_x(0,t) &= &A=w_x(0,t)+v_x(0,t)=A+v_x(0,t) &\Rightarrow &v_x(0,t)=0\\
u_x(L,t) &= &B=w_x(L,t)+v_x(L,t)=B+v_x(L,t) &\Rightarrow &v_x(L,t)=0\\
u(x,0) &= &g(x)=w(x,0)+v(x,0) &\Rightarrow
&v(x,0)=g(x)-w(x,0)\end{array}. \label{eg171veq}
\end{eqnarray}
$$(eg171veq)

Equations {eq}`eg171veq` represent the homogeneous Neumann BVP
seen previously. Therefore

$$
\begin{eqnarray}
u(x,t) & = & \frac{(B-A)}{2L} x^2+Ax+\alpha
^2\left(\frac{B-A}{L}\right)
t+\frac{a_0}{2}\\
& &\hspace{1in}+\sum\limits_{n=1}^\infty a_n\cos\left(\frac{n\pi
x}{L}\right)\{\rm\ e\}^{-\alpha^2\left(\frac{n\pi x}{L}\right) t}\nonumber
\end{eqnarray}
$$(ref-heat-inhomogeneous-42)

where

$$
\begin{eqnarray}
a_n=\frac{2}{L}\int\limits_0^L\left\{
g(x)-\left[\frac{(B-A)}{2L}x^2+Ax\right]\right\}\cos\left(\frac{n\pi
x}{L}\right)\, dx.
\end{eqnarray}
$$(ref-heat-inhomogeneous-43)
````
