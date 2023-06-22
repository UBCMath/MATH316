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

# Heat Conduction With Time Dependent Boundary Conditions Using Eigenfunction Expansions

The ultimate goal of this lecture is to demonstrate a method to
solve heat conduction problems in which there are time dependent
boundary conditions. The idea is to construct the simplest possible
function, $w(x,t)$ say, that satisfies the inhomogeneous,
time-dependent boundary conditions. The solution $u(x,t)$ that we
seek is then decomposed into a sum of $w(x,t)$ and another function
$v(x,t)$, which satisfies the homogeneous boundary conditions. When
these two functions are substituted into the heat equation, it is
found that $v(x,t)$ must satisfy the heat equation subject to a
source that can be time dependent. As in Lecture 19, this forced
heat conduction equation is solved by the method of eigenfunction
expansions.

```{admonition} Key Concepts
Time-dependent Boundary conditions,
distributed sources/sinks, Method of Eigenfunction Expansions.
```

## Heat Conduction Problems with Time Dependent Boundary Conditions

### Inhomogeneous Derivative Boundary Conditions Using Eigenfunction Expansions

````{prf:example}
:label: example-heat-eigenfunctions-0

Let us revisit the problem with inhomogeneous derivative BC {prf:ref}`example-heat-inhomogeneous-3`
- but we will now use Eigenfunction Expansions.

$$
\begin{eqnarray}
u_t & = & \alpha ^2u_{xx}\qquad 0<x<L,\quad t>0\\
\mbox{BC: } u_x(0,t) & = & A\quad u_x(L,t)=B\\
\mbox{IC: } u(x,0) &=& g(x)
\end{eqnarray}
$$(ref-heat-eigenfunctions-0)

First look for a function of the form $h(x)=ax^2+bx$ that satisfies
the inhomogeneous BC:

$$
\begin{eqnarray*}
h(x) &=& ax^2+bx,\quad h_x(x)=2ax+b\\
h_x(0) &=& b=A\quad h_x(L)=2aL+A=B\Rightarrow a=(B-A)/2L\\
h(x) &=& \left(\frac{B-A}{2L}\right) x^2+Ax.
\end{eqnarray*}
$$(ref-heat-eigenfunctions-1)

Now let

$$
\begin{eqnarray*}
u(x,t)=h(x)+v(x,t).
\end{eqnarray*}
$$(ref-heat-eigenfunctions-2)

Substitute into the PDE:

$$
\begin{eqnarray*}
u_t={\big( h(x)\!\!\!\!\!\!\!\!\!\nearrow
+v(x,t)\big)}_t=\alpha^2u_{xx}=\alpha^2 {\big(
h(x)+v(x,t)\big)}_{xx}=\alpha^2\cdot 2a+\alpha^2 v_{xx}.
\end{eqnarray*}
$$(ref-heat-eigenfunctions-3)

Therefore

$$
\begin{eqnarray}
v_t=\alpha^2 v_{xx}+2a\alpha^2 \label{lec21ee1eq1}
\end{eqnarray}
$$(lec21ee1eq1)

$$
\begin{eqnarray}
A=u_x(0,t)=h_x(0)+v_x(0,t)=A+V_x(0,t) &\Rightarrow & v_x(0,t)=0\\
B=u_x(L,t)=h_x(L)+V_x(L,t)=B+V_x(L,t) &\Rightarrow & v_x(L,t)=0\\
g(x)=u(x,0)=h(x)+v(x,0) &\Rightarrow &
v(x,0)=g(x)-h(x).\phantom{\Rightarrow\Rightarrow}
\label{lec21ee1eq2}
\end{eqnarray}
$$(lec21ee1eq2)

We now use an Eigenfunction Expansion to solve the BVP
{eq}`lec21ee1eq1`-{eq}`lec21ee1eq2`. Because of the homogeneous
Neumann BC we assume an expansion of the form

$$
\begin{eqnarray*}
v(x,t) &=& \hat{v}_0(t)/2+\sum\limits_{n=1}^\infty \hat{v}_n(t)\cos\left(\frac{n\pi x}{L}\right)\\
v_t &=& \dot{\hat{v}}_0  (t)/2+\sum\limits_{n=1}^\infty \dot{\hat{v}}_n (t)\cos\left(\frac{n\pi x}{L}\right)\\
v_x &=& \sum\limits_{n=1}^\infty \hat{v}_n(t)\left\{
-\left(\frac{n\pi}{L}\right)\right\}\sin\left(\frac{n\pi
x}{L}\right) ,v_{xx}=\sum\limits_{n=1}^\infty \hat{v}_n(t)\left\{
-\left(\frac{n\pi}{L}\right)^2\right\}\cos\left(\frac{n\pi
x}{L}\right) .
\end{eqnarray*}
$$(ref-heat-eigenfunctions-6)

We also expand the inhomogeneous term in {eq}`lec21ee1eq1` in terms of the
Eigenfunctions:

$$
\begin{eqnarray*}
2a\alpha^2 =a_0/2+\sum\limits_{n=1}^\infty a_n\cos\left(\frac{n\pi
x}{L}\right)\quad a_0=4a\alpha^2, a_n=0\quad n\geq 1.
\end{eqnarray*}
$$(ref-heat-eigenfunctions-7)

Therefore

$$
\begin{eqnarray*}
0=v_t-\alpha^2 v_{xx}-2a\alpha^2 =\dot{\hat{v}}_0 (t)/2-2a\alpha^2
+\sum\limits_{n=1}^\infty \left\{ \dot{\hat{v}}_n +\alpha^2
\left(\frac{n\pi}{L}\right)^2 \hat{v}_n\right\} \cos\left(\frac{n\pi
x}{L}\right) .
\end{eqnarray*}
$$(ref-heat-eigenfunctions-8)

Therefore

$$
\begin{eqnarray*}
\dot{\hat{v}}_0(t)&=& 4a\alpha^2 \Rightarrow\hat{v}_0(t)=4a\alpha^2 t+c_0\\
\dot{\hat{v}}_n(t) &=& -\alpha^2\left(\frac{n\pi}{L}\right)^2
\hat{v}_n\Rightarrow \hat{v}_n(t)=\hat{v}_n(0)
\{\rm\ e\}^{-\alpha^2\left(\frac{n\pi}{L}\right)^2 t}.
\end{eqnarray*}
$$(ref-heat-eigenfunctions-9)

Therefore

$$
\begin{eqnarray*}
v(x,t)=\frac{4a\alpha^2t+c_0}{2}+\sum\limits_{n=1}^\infty \hat{v}_n (0)
\{\rm\ e\}^{-\alpha^2 \left(\frac{n\pi}{L}\right)^2t}\cos\left(\frac{n\pi
x}{L}\right) .
\end{eqnarray*}
$$(ref-heat-eigenfunctions-10)

$$
\begin{eqnarray*}
g(x)-h(x) &=& g(x)-\left\{\left( \frac{B-A}{2L}\right)
x^2+Ax\right\}
=v(x,0)=\frac{c_0}{2}+\sum\limits_{n=1}^\infty \hat{v}_n(0)\cos\left(\frac{n\pi x}{L}\right)\\
c_0 &=&\frac{2}{L}\int\limits_0^1 \left[ g(x)-\left\{\left(\frac{B-A}{2L}\right) x^2+Ax\right\}\right]\, dx\\
\hat{v}_n(0) &=& \frac{2}{L}\int\limits_0^1 \left[
g(x)-\left\{\left(\frac{B-A}{2L}\right) x^2+Ax\right\} \right]
\cos\left(\frac{n\pi x}{L}\right)\, dx.
\end{eqnarray*}
$$(ref-heat-eigenfunctions-11)

Thus

$$
\begin{eqnarray*}
u(x,t)=\left(\frac{B-A}{2L}\right) x^2+Ax+2a\alpha^2
t+\frac{c_0}{2}+\sum\limits_{n=1}^\infty \hat{v}_n
(0)\{\rm\ e\}^{-\alpha^2\left(\frac{n\pi}{L}\right)^2 t}\cos\left(\frac{n\pi
x}{L}\right)
\end{eqnarray*}
$$(ref-heat-eigenfunctions-12)

which is identical to the solution obtained in {prf:ref}`example-heat-inhomogeneous-3`.

````

### Time-Dependent Boundary Conditions Using Eigenfunction Expansions

````{prf:example} Time Dependent Boundary Conditions - General Case
:label: example-heat-eigenfunctions-1 

$$
\begin{eqnarray}
u_t & = & \alpha^2 u_{xx},\qquad 0<x<L\nonumber\\
\mbox{BC: } u(0,t) & = & \phi_0 (t)\quad u(L,t)=\phi_1 (t)\quad
\nonumber\\
\mbox{IC: } u(x,0) & = & f(x).
\end{eqnarray}
$$(ref-heat-eigenfunctions-13)

```{figure} ../img/heat/time_dependent_dirichlet.png
:name: time_dependent_dirichlet
:alt: time_dependent_dirichlet
:align: center

Bar subject to a time-dependent Dirichlet BC.
```

Let \framebox{$\displaystyle w(x,t)=\phi_0
(t)+x\left(\frac{\phi_1(t)-\phi_0(t)}{L}\right)$} $\Rightarrow
w(0,t)=\phi_0 (t)$; $w(L,t)=\phi_1(t)$. Now let
$u(x,t)=w(x,t)+v(x,t)$. Then

$$
\begin{eqnarray}
w_t+v_t & = & \alpha^2 (w\!\!\!\!\!\nearrow_{\!\!\!\! xx}+v_{xx})\nonumber\\
v_t & = & \alpha^2 v_{xx}-w_t\qquad w_t=\dot{\phi}_0
+\frac{x}{L}(\dot{\phi}_1-\dot{\phi}_0)\nonumber\\
\mbox{BC: } u(0,t) & =& \phi_0(t)=w(0,t)+v(0,t)=\phi_0
(t)+v(0,t)\Rightarrow v(0,t)=0\nonumber\\
u(L,t) & = & \phi_1(t)=w(L,t)+v(L,t)=\phi_1(t)+v(L,t)\Rightarrow
v(L,t)=0\nonumber\\
\mbox{IC: }u(x,0) & = & f(x)=w(x,0)+v(x,0)\Rightarrow
v(x,0)=f(x)-w(x,0).
\end{eqnarray}
$$(ref-heat-eigenfunctions-14)

Thus we need to solve the following BVP for $v(x,t)$:

$$
\begin{eqnarray}
v_t & = & \alpha^2 v_{xx}-w_t\nonumber\\
\mbox{BC: }v(0,t) & = & 0\quad v(L,t)=0\\
\mbox{IC: }v(x,0) & = & f(x)-w(x,0).\nonumber
\end{eqnarray}
$$(ref-heat-eigenfunctions-15)

Now $v(x,t)$ can be found using an eigenfunction expansion. The
eigenfunctions and eigenvalues associated with the Dirichlet B C are

$$
\begin{eqnarray*}
\lambda_n &=&\left(\frac{n\pi}{L}\right)\quad n=1,2,\ldots\quad X_n(x)=\sin (\lambda_n x)\\
\mbox{let }S(x,t)&=& -w_t \\
&=& - {(\dot{\phi}_1
-\dot{\phi}_0)\left(\frac{x}{L}\right) -\dot{\phi}_0} \\
&=& \sum\limits_{n=1}^\infty \hat{S}_n (t)\sin (\lambda_n x)\\
\mbox{and }v(x,t)&=&\sum\limits_{n=1}^\infty \hat{v}_n(t)\sin (\lambda_n x)\\
\mbox{then }v_t&=&\sum\limits_{n=1}^\infty \dot{\hat{v}}_n (t)\sin
(\lambda_nx)\mbox{ and }
  v_{xx}=\sum\limits_{n=1}^\infty \hat{v}_n (t)\left\{ -\lambda_n^2\right\}\sin (\lambda_n x)\\
\mbox{thus }0&=&v_t -\alpha^2 v_{xx}-S(x,t)
\end{eqnarray*}
$$(ref-heat-eigenfunctions-16)

Therefore

$$
\begin{eqnarray}
0=\sum\limits_{n=1}^\infty\left\{
\dot{\hat{v}}_n+\alpha^2\lambda_n^2
\hat{v}_n-\hat{S}_n(t)\right\}\sin (\lambda_n x) \label{eigeneq}
\end{eqnarray}
$$(eigeneq)

Since the eigenfunctions are linearly independent it follows that
$\{\,\,\} =0$ in {eq}`eigeneq` or

$$
\begin{eqnarray}
\frac{d\hat{v}_n}{dt} +\alpha^2 \lambda_n^2\hat{v}_n=\hat{S}_n (t)
\label{eigenode}
\end{eqnarray}
$$(eigenode)

but {eq}`eigenode` is just a first order linear ODE with an
integrating factor

$$
\begin{eqnarray*}
F(t)=e^{\alpha^2\lambda_n^2t}
\end{eqnarray*}
$$(ref-heat-eigenfunctions-19)

Thus

$$
\begin{eqnarray*}
\frac{d}{dt}\left( e^{\alpha^2\lambda_n^2t}\hat{v}_n(t)\right)
=e^{\alpha^2\lambda_n^2t}\hat{S}_n(t)
\end{eqnarray*}
$$(ref-heat-eigenfunctions-20)

Integrating we obtain 

$$
\begin{eqnarray*}
e^{\alpha^2\lambda_n^2t}\hat{v}_n(t)=\int\limits_0^t
e^{\alpha^2\lambda_n^2\tau}\hat{S}_n(\tau )\, d\tau +c_n
\end{eqnarray*}
$$(ref-heat-eigenfunctions-21)

or

$$
\begin{eqnarray*}
\hat{v}_n(t)=\int\limits_0^t e^{-\alpha^2\lambda_n^2 (t-\tau
)}\hat{S}_n(\tau )\, d\tau +e^{-\alpha\lambda_n^2 t}c_n
\end{eqnarray*}
$$(ref-heat-eigenfunctions-22)

Thus

$$
\begin{eqnarray*}
v(x,t)=\sum\limits_{n=1}^\infty \left\{\int\limits_0^t
e^{-\alpha^2\lambda_n^2 (t-\tau )}\hat{S}_n(\tau )\, d\tau
+e^{-\alpha^2\lambda_n^2t}c_n\right\}\sin (\lambda_n x)
\end{eqnarray*}
$$(ref-heat-eigenfunctions-23)

All we need to do to complete the solution of this problem is to
determine the coefficients $c_n$. These we obtain from the initial
condition as follows

$$
\begin{eqnarray*}
g(x)-\left[\left\{\phi_1 (0)-\phi_0
(0)\right\}\left(\frac{x}{L}\right) +\phi_0 (0)\right]
=\sum\limits_{n=1}^\infty c_n\sin (\lambda_n x)
\end{eqnarray*}
$$(ref-heat-eigenfunctions-24)

But this is just a Fourier sine series in which

$$
\begin{eqnarray*}
c_n=\frac{2}{L}\int\limits_0^L\left( g(x)-\left[\left\{\phi_1
(0)-\phi_0 (0)\right\}\left(\frac{x}{L}\right) -\phi_0
(0)\right]\right)\sin\left(\frac{n\pi x}{L}\right)\, dx
\end{eqnarray*}
$$(ref-heat-eigenfunctions-25)

Finally

$$
\begin{eqnarray*}
u(x,t)=\left(\phi_1(t)-\phi_0 (t)\right)\left(\frac{x}{L}\right)
+\phi_0 (t)+\sum\limits_{n=1}^\infty\left\{\int\limits_0^t
e^{-\alpha^2\lambda_n^2 (t-\tau )}\hat{S}_n(\tau )\, d\tau
+e^{-\alpha^2\lambda_n^2t}c_n\right\}\sin\lambda_n x.
\end{eqnarray*}
$$(ref-heat-eigenfunctions-26)

 __Specific case:__ 
 
Let $\phi_0(t)=At$, $\phi_1(t)=0$, and $f(x)=0$. In this case

$$
\begin{eqnarray}
w(x,t)=At+\frac{x}{L}(0-At)=At\left( 1-\frac{x}{L}\right) .
\end{eqnarray}
$$(ref-heat-eigenfunctions-27)

$$
\begin{eqnarray}
u_t & = & \alpha^2 u_{xx}\qquad 0<x<L\nonumber\\
\mbox{BC: }u(0,t) & = & At\quad u(L,t)=0\\
\mbox{IC: }u(x,t) & = & 0.\nonumber
\end{eqnarray}
$$(ref-heat-eigenfunctions-28)

Let $u(x,t)=w(x,t)+v(x,t)$ where $\displaystyle w(x,t)=At\left(
1-\frac{x}{L}\right)$. Then

$$
\begin{eqnarray}
v_t & = & \alpha^2 v_{xx}-A\left( 1-\frac{x}{L}\right)\nonumber\\
v(0,t) & = & 0=v(L,t)\\
v(x,0) & = & 0.\nonumber
\end{eqnarray}
$$(ref-heat-eigenfunctions-29)

Let

$$
\begin{eqnarray}
s(x,t ) & = & -A\left( 1-\frac{x}{L}\right)
=\sum\limits_{n=1}^\infty \hat{s}_n(t)\sin\left(\frac{n\pi
x}{L}\right)\nonumber\\
\hat{s}_n & = & \frac{2}{L}\int\limits_0^L
A\left(\frac{x}{L}-1\right)\sin\left(\frac{n\pi x}{L}\right)\, dx\\
& = & -\frac{2A}{n\pi}.\nonumber
\end{eqnarray}
$$(ref-heat-eigenfunctions-30)

Now let

$$
\begin{eqnarray}
v(x,t) & = & \sum\limits_{n=1}^\infty
\hat{v}_n(t)\sin\left(\frac{n\pi x}{L}\right)\nonumber\\
\\
v_t & = & \sum\limits_{n=1}^\infty
\dot{\hat{v}}_n(t)\sin\left(\frac{n\pi x}{L}\right), \quad
v_{xx}=-\sum\limits_{n=1}^\infty \hat{v}_n(t){\left(\frac{n\pi
}{L}\right)}^2\sin\left(\frac{n\pi x}{L}\right) .\nonumber
\end{eqnarray}
$$(ref-heat-eigenfunctions-31)

Therefore

$$
\begin{eqnarray}
0=v_t-\alpha^2 v_{xx}-s(x,t)=\sum\limits_{n=1}^\infty\left\{
  \dot{\hat{v}}_n(t)+\alpha^2 {\left(\frac{n\pi}{L}\right)}^2
  \hat{v}_n+\frac{2A}{n\pi}\right\}\sin\left(\frac{
  n\pi x}{L}\right) .
\end{eqnarray}
$$(ref-heat-eigenfunctions-32)

Therefore

$$
\begin{eqnarray}
\dot{\hat{v}}_n(t)+\alpha^2 {\left(\frac{n\pi}{L}\right)}^2
\hat{v}_n(t) & = & -\frac{2A}{n\pi}\\
\left( \{\rm\ e\}^{+\alpha^2
{\left(\frac{n\pi}{L}\right)}^2t}\hat{v}_n(t)\right) & = &
-\frac{2A}{n\pi}\{\rm\ e\}^{\alpha^2 {\left(\frac{n\pi}{L}\right)}^2 t}\\
\{\rm\ e\}^{\alpha^2\left(\frac{n\pi}{L}\right)^2 t}\hat{v}_n(t) & = &
-\frac{2AL^2}{\alpha^2(n\pi )^3}
  \{\rm\ e\}^{\alpha^2\left(\frac{n\pi}{L}\right)^2 t}+B_n\\
\hat{v}_n(t) & = &
  -\frac{2AL^2}{\alpha^2 (n\pi )^3}+B_n
    \{\rm\ e\}^{-\alpha^2\left(\frac{n\pi}{L}\right)^2 t}\\
0=\hat{v}_n(0) & = & -\frac{2AL^2}{\alpha^2(n\pi )^3}+B_n.
\end{eqnarray}
$$(ref-heat-eigenfunctions-33)

Therefore

$$
\begin{eqnarray}
\hat{v}_n(t) & = & \frac{2AL^2}{\alpha^2(n\pi)^3}
   \left(\{\rm\ e\}^{-\alpha^2\left(\frac{n\pi}{L}\right)^2 t}-1 \right) .
\end{eqnarray}
$$(ref-heat-eigenfunctions-34)

Therefore

$$
\begin{eqnarray}
u(x,t)=At\left(1-\frac{x}{L}\right) +
\frac{2AL^2}{\pi^3\alpha^2}\sum\limits_{n=1}^\infty
   \frac{(\{\rm\ e\}^{-\alpha^2\left(\frac{n\pi}{L}\right)^2 t}-1)}{n^3}
      \sin\left(\frac{n\pi x}{L}\right) .
\end{eqnarray}
$$(ref-heat-eigenfunctions-35)
````

### Summary of Guesses for $w(x,t)$ to Remove Different Inhomogeneous Boundary Conditions

Consider the following heat equation subject to a loss represented by $-\gamma u$ and a
source $S(x,t)$:

$$
\begin{eqnarray*}
u_t=\alpha^2 u_{xx}-\gamma u+S(x,t)
\end{eqnarray*}
$$(ref-heat-eigenfunctions-36)

#### Mixed BC I

$$
\begin{eqnarray*}
u(0,t)=\phi_0(t),\quad u_x(L,t)=\phi_1(t),\quad w=\phi_0+\phi_1x
\end{eqnarray*}
$$(ref-heat-eigenfunctions-37)

#### Mixed BC II

$$
\begin{eqnarray*}
u_x(0,t)=\phi_0(t),\quad u(L,t)=\phi_1(t),\quad w=(\phi_1
-\phi_0L)+\phi_0x
\end{eqnarray*}
$$(ref-heat-eigenfunctions-38)

#### Dirichlet BC

$$
\begin{eqnarray*}
u(0,t)&=&\phi_0(t),\quad u(L,t)=\phi_1(t),\quad
w=\phi_0+(\phi_1-\phi_0)\frac{x}{L}
\end{eqnarray*}
$$(ref-heat-eigenfunctions-39)

#### Neumann BC

$$
\begin{eqnarray*}
u_x(0,t)=\phi_0(t),\quad u_x(L,t)=\phi_1(t),\quad
w=\phi_0x+(\phi_1-\phi_0)\frac{x^2}{2L}
\end{eqnarray*}
$$(ref-heat-eigenfunctions-40)

Let $u(x,t)=w(x,t)+v(x,t)$

$$
\begin{eqnarray*}
w_t+v_t&=&\alpha^2(w_{xx}+v_{xx})-\gamma (w+v)+S(x,t)\\
v_t&=&\alpha^2 v_{xx}-\gamma v+\left\{\alpha^2w_{xx}-\gamma
w-w_t\right\} +S(x,t)
\end{eqnarray*}
$$(ref-heat-eigenfunctions-41)
