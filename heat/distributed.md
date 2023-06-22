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

# Heat Conduction With Distributed Sources/Sinks

In this lecture we consider heat conduction problems in which there
is a distributed source or sink function $s(x,t)$ that applies
throughout the domain. We first consider the case in which the
source is time-independent, i.e., $s(x,t)=s(x)$. In this case the
effect of the source can be dealt with entirely by determining an
appropriate steady-state solution. Using this particular solution,
we can reduce the problem to one of the standard homogeneous
boundary value problems we encountered when we introduced separation
of variables. Secondly, we consider a fully time-dependent source.
In this case we have to resort to the method of eigenfunction
expansions.

```{admonition} Key Concepts
Distributed sources or sinks,
Particular Solutions, Steady state Solutions; Separation of
variables, Eigenvalues and Eigenfunctions, Method of Eigenfunction
Expansions.
```

## Heat Conduction Problems With Distributed Sources

### Heat Conduction Problems With Distributed Time-Independent Sources

````{prf:example} A bar with an external heat source $s(x)=x$.
:label: example-heat-distributed-0 

$$
\begin{eqnarray}
u_t & = & \alpha^2 u_{xx}+x\qquad 0<x<L\\
\mbox{BC: }u(0,t) & = & 0\quad u(L,t)=B\\
\mbox{IC: }u(x,0) & = & g(x).
\end{eqnarray}
$$(ref-heat-distributed-0)

```{figure} ../img/heat/time_independent_sources.png
:name: time_independent_sources.png
:alt: time_independent sources.png
:align: center

Bar subject to a time-independent heat source distributed
along its length with inhomogeneous Dirichlet BC
```

Steady state problem $u_t=0$:

$$
\begin{eqnarray}\begin{array}{l}
0=\alpha^2 u_\infty^{\prime\prime}+x\\
u_\infty (0)=0\quad u_\infty (L)=B\end{array}
\end{eqnarray}
$$(ref-heat-distributed-1)

$$
\begin{eqnarray}\begin{array}{l}
\displaystyle u_\infty^{\prime\prime}=-\frac{x}{\alpha^2}\quad
u_\infty^\prime
  =-\frac{x^2}{2\alpha^2}+a\quad u_\infty =-\frac{x^3}{6\alpha^2}+ax+b\\
\displaystyle u_\infty (0)=b=0\quad u_\infty
(L)=-\frac{L^3}{6\alpha^2}+aL=B\Rightarrow
a=\frac{B}{L}+\frac{L^2}{6\alpha^2}\end{array}
\end{eqnarray}
$$(ref-heat-distributed-2)

Therefore

$$
\begin{eqnarray}
u_\infty
(x)=-\frac{x^3}{6\alpha^2}+\left(\frac{B}{L}+\frac{L^2}{6\alpha^2}\right)
x=x\left\{\frac{B}{L}+\frac{1}{6\alpha^2}(L^2-x^2)\right\} .
\end{eqnarray}
$$(ref-heat-distributed-3)

Let $u(x,t)=u_\infty (x)+v(x,t)$.

$$
\begin{eqnarray}\begin{array}{lclcl}
u_t=\alpha^2 u_{xx}+x &\Rightarrow &(u_\infty\!\!\!\!\nearrow +v)_t
  =\alpha^2(u\!\!\!\downarrow _\infty +v)_{xx}+x\!\!\!\!\!\downarrow &\Rightarrow &v_t=\alpha^2 v_{xx}\\
u(0,t)=0 &\Rightarrow &u_\infty \!\!\!\!\nearrow (0)+v(0,t)=0 &\Rightarrow &v(0,t)=0\\
u(L,t)=B &\Rightarrow &u_\infty (L)+v(L,t)=B &\Rightarrow &v(L,t)=0\\
u(x,0)=g(x) &\Rightarrow &u_\infty (x)+v(x,0)=g(x) &\Rightarrow
&v(x,0)=g(x)-u_\infty (x).\end{array}
\end{eqnarray}
$$(ref-heat-distributed-4)

Separation of variables yields:

$$
\begin{eqnarray}
v(x,t)=\sum\limits_{n=1}^\infty
b_n\{\rm\ e\}^{-{\left(\frac{n\pi}{L}\right)}^2\alpha^2 t}
   \sin\left(\frac{n\pi x}{L}\right)
\end{eqnarray}
$$(ref-heat-distributed-5)

where

$$
\begin{eqnarray}
b_n=\frac{2}{L}\int\limits_0^L\left\{ g(x)-u_\infty (x)\right\}
   \sin\left(\frac{n\pi x}{L}\right)\, dx.
\end{eqnarray}
$$(ref-heat-distributed-6)

Therefore

$$
\begin{eqnarray}
u(x,t) & = & x\left\{\frac{B}{L}+\frac{1}{6\alpha^2}
(L^2-x^2)\right\}
   +\sum\limits_{n=1}^\infty b_n\{\rm\ e\}^{-\alpha^2 {\left(\frac{n\pi}{L}\right)}^2 t}
   \sin\left(\frac{n\pi x}{L}\right)\nonumber\\
& &\hspace{.75in}\uparrow\hspace{1.5in}\uparrow\\
& &\hspace{.75in}\mbox{steady}\hspace{1in}\mbox{transient}\nonumber
\end{eqnarray}
$$(ref-heat-distributed-7)

```{note}
$$
\begin{eqnarray}
\lim\limits_{x\rightarrow\infty}
u(x,t)=x\left\{\frac{B}{L}+\frac{1}{6\alpha^2} (L^2-x^2)\right\} .
\end{eqnarray}
$$(ref-heat-distributed-8)
```
````

### Distributed time-dependent heat sources/sinks - eigenfunction expansions

````{prf:example} A Bar with a Time-Varying External Heat Source
:label: example-heat-distributed-1 

$$
\begin{eqnarray}
u_t & = & \alpha^2 u_{xx}+\{\rm\ e\}^{-t}\sin\left(\frac{2\pi x}{L}\right)\quad 0<x<L,\quad t>0\\
\mbox{BC: }u(0,t) & = & 0;\quad u(L,t)=L\\
\mbox{IC: }u(x,0) & = & x.
\end{eqnarray}
$$(ref-heat-distributed-9)

Consider the function $w(x)=x$ which satisfies the BC as well as the
homogeneous version of the PDE.

Now let $u(x,t)=w(x)+v(x,t)$.

$$
\begin{eqnarray}
u_t & = & (w\!\!\!\!\!\!\nearrow +v)_t =
\alpha^2(w\!\!\!\!\!\!\nearrow +v)_{xx}
  +\{\rm\ e\}^{-t}\sin\left(\frac{2\pi x}{L}\right)\\
   & &\qquad\Rightarrow v_t
   =\alpha^2 v_{xx}+\{\rm\ e\}^{-t}\sin\left(\frac{2\pi x}{L}\right)\\
u(0,t) & = & w(0)\!\!\!\!\!\!\!\!\!\!\!\nearrow +v(0,t)=0\Rightarrow v(0,t)=0\\
u(L,t) & = & w(L)\!\!\!\!\!\!\!\!\!\!\!\nearrow +v(L,t)=L\!\!\!\!\!\!\!\nearrow \Rightarrow v(L,t)=0\\
x=u(x,0) & = & w(x)+v(x,0)=x+v(x,0)\Rightarrow v(x,0)=0.
\end{eqnarray}
$$(ref-heat-distributed-10)

Now assume that $\displaystyle v(x,t)=\sum\limits_{n=1}^\infty
\hat{v}_n(t)\sin\left(\frac{n\pi x}{L}\right)$.

$$
\begin{eqnarray}
\frac{\partial v}{\partial t}=\sum\limits_{n=1}^\infty
\frac{d\hat{v}_n}{dt}(t)\sin\left(\frac{n\pi x}{L}\right),\; \;
\frac{\partial^2 v}{\partial x^2}=\sum\limits_{n=1}^\infty \hat{v}_n
(t)\left\{ -{\left(\frac{n\pi}{L}\right) }^2 \sin\left(\frac{n\pi
x}{L}\right) \right\} .
\end{eqnarray}
$$(ref-heat-distributed-11)

Therefore

$$
\begin{eqnarray}
v_t-\alpha^2 v_{xx}-\{\rm\ e\}^{-t}\sin\left(\frac{2\pi x}{L}\right)=\sum\limits_{n=1}^\infty
\left\{\frac{d\hat{v}_n}{dt}+\alpha^2
{\left(\frac{n\pi}{L}\right)}^2\hat{v}_n-\{\rm\ e\}^{-t}\delta_{2n}\right\}\sin\left(\frac{n\pi
x}{L}\right) =0.
\end{eqnarray}
$$(ref-heat-distributed-12)

Therefore

$$
\begin{eqnarray}
\frac{d\hat{v}_n}{dt} +\alpha^2 {\left(\frac{n\pi}{L}\right)}^2
\hat{v}_n & = &\{\rm\ e\}^{-t}\delta_{2n}\\
\frac{d}{dt}\left[ \{\rm\ e\}^{\alpha^2 {\left(\frac{n\pi}{L}\right)}^2 t}
\hat{v}_n\right]
  & = & \{\rm\ e\}^{\left[\alpha^2 {\left(\frac{n\pi}{L}\right)}^2-1\right] t}\delta_{2n}.
\end{eqnarray}
$$(ref-heat-distributed-13)

Therefore

$$
\begin{eqnarray}
\{\rm\ e\}^{-\alpha^2 \left( \frac{n\pi}{L}\right)^2 t}\hat{v}_n
  & = & \frac{\{\rm\ e\}^{\left[\alpha^2 {\left(\frac{n\pi}{L}\right)}^2 -1\right] t}}
  {\alpha^2 {\left(\frac{n\pi}{L}\right)}^2 -1}
   \delta_{2n}+c_n\qquad c_n\mbox{ arbitrary}\\
\hat{v}_n(t) & = & \frac{\{\rm\ e\}^{-t}\delta_{2n}}
{\alpha^2{\left(\frac{n\pi}{L}\right)}^2 -1}
   +\{\rm\ e\}^{-\alpha^2
{\left(\frac{n\pi}{L}\right)}^2 t}c_n\\
v(x,0) & = & 0\Rightarrow \hat{v}_n(0)=0=\frac{\delta_{2n}}
   {\alpha^2 {\left(\frac{n\pi}{L}\right)}^2 -1}
   +c_n\Rightarrow\\
     & &\hspace{1.5in} c_n=\left\{ \begin{array}{ll}
       0 &n\ne 2\\ -\frac{1}{\alpha^2 {\left(\frac{2\pi}{L}\right)}^2
       -1} &n=2\end{array}\right.\nonumber\\
v(x,t) & = & \frac{1}{\alpha^2{\left(\frac{2\pi}{L}\right)}^2 -1}
   \left\{
     \{\rm\ e\}^{-t}-\{\rm\ e\}^{-\alpha^2 {\left(\frac{2\pi}{L}
       \right)}^2 t}\right\}
   \sin\left(\frac{2\pi x}{L}\right)\\
u(x,t) & = & x+v(x,t)=x+ \left(
   \frac{\{\rm\ e\}^{-t} -
   \{\rm\ e\}^{-\alpha^2\left(\frac{2\pi}{L}\right)^2 t}}
{\alpha^2 \left(\frac{2\pi}{L}\right)^2-1}
    \right)
    \sin\left(\frac{2\pi x}{L}\right) . \nonumber
\end{eqnarray}
$$(ref-heat-distributed-14)
````

````{prf:example} A bar with a general external heat source $s(x,t)$
:label: example-heat-distributed-2 

$$
u_t = \alpha^2 u_{xx}+s(x,t) \label{eg172a}
$$(eg172a)

$$
\mbox{BC: } u(0,t) = A\quad u(L,t)=B \label{eg172b}
$$(eg172b)

$$
\mbox{IC: }u(x,0) = f(x). \label{eg172c}
$$(eg172c)

```{figure} ../img/heat/general_heat_source.png
:name: general_heat_source
:alt: general_heat_source
:align: center

Bar subject to a time dependent heat source distributed
along its length with inhomogeneous Dirichlet BC
```

We look for a particular solution: $w(x,t)$ by expanding $s(x,t)$ as
a Sine Series. Note that the sine functions are the eigenfunctions
that correspond to the homogeneous form of the BC in {eq}`eg172b`.
Thus if we add $w(x,t)$ to a solution of
{eq}`eg172a`-{eq}`eg172b` without the source (i.e. with
$s(x,t)=0$) we will not affect the BC.

__Eigenfunction Expansion:__

Let

$$
\begin{equation}
s(x,t)=\sum\limits_{n=1}^\infty \hat{s}_n(t)\sin\left(\frac{n\pi
x}{L}\right)
\end{equation}
$$(ref-heat-distributed-16)

where

$$
\begin{equation}
\hat{s}_n(t)=\frac{2}{L}\int\limits_0^L s(x,t)\sin\left(\frac{n\pi
x}{L}\right)\, dx.
\end{equation}
$$(ref-heat-distributed-17)

If we assume

$$
\begin{equation}
w(x,t)=\sum\limits_{n=1}^\infty \hat{w}_n(t)\sin\left(\frac{n\pi
x}{L}\right)
\end{equation}
$$(ref-heat-distributed-18)

then

$$
\begin{eqnarray}
w_t & = & \sum\limits_{n=1}^\infty \hat{w}_n^\prime (t)
  \sin\left(\frac{n\pi x}{L}\right)\\
w_{xx} & = & -\sum\limits_{n=1}^\infty \hat{w}_n
  {\left(\frac{n\pi}{L}\right)}^2 \sin\left(\frac{n\pi x}{L}\right) .
\end{eqnarray}
$$(ref-heat-distributed-19)

Therefore substituting these expansions into $w_t=\alpha^2
w_{xx}+s(x,t)$ we obtain:

$$
\begin{equation}
\sum\limits_{n=1}^\infty \left\{ \hat{w}_n^\prime +\alpha^2 \left(
\frac{n\pi}{L}\right)^2 \hat{w}_n- \hat{s}_n (t)\right\}
\sin\left(\frac{n\pi x}{L}\right) =0.
\end{equation}
$$(ref-heat-distributed-20)

Therefore

$$
\begin{equation}
\hat{w}_n^\prime (t)= -\alpha^2 {\left(\frac{n\pi}{L}\right)}^2
\hat{w}_n (t) +\hat{s}_n(t).
\end{equation}
$$(ref-heat-distributed-21)

This is a linear 1st order ODE with integrating factor $
\{\rm\ e\}^{\alpha^2 {\left(\displaystyle\frac{n\pi}{L}\right)}^2 t}$.

Therefore

$$
\begin{equation}
w_n(t)=\int\limits_0^t \{\rm\ e\}^{-
\alpha^2{\left(\frac{n\pi}{L}\right)}^2 (t-\tau )}
   \hat{s}_n(\tau )\, d\tau +c_n\{\rm\ e\}^{-\alpha^2 {\left(\frac{n\pi}{L}\right)}^2 t}
\end{equation}
$$(ref-heat-distributed-22)

where the $c_n$ are arbitrary constants. Since we are only looking
for a particular solution we choose $c_n\equiv 0$.

Therefore

$$
\begin{eqnarray}
w(x,t)=\sum\limits_{n=1}^\infty \left( \int\limits_0^t
  \{\rm\ e\}^{-\alpha^2 \left(\frac{n\pi}{L}\right)^2 (t-\tau )}
  \hat{s}_n (\tau )
\, d\tau \right) \sin\left(\frac{n\pi x}{L}\right) .
\end{eqnarray}
$$(ref-heat-distributed-23)

__Exploiting Linearity via Superposition:__

Now that we have a particular solution we exploit the fact that
the Problem {eq}`eg172a`-{eq}`eg172b` is linear and use
superposition.

Let

$$
\begin{eqnarray}
u(x,t) & = & w(x,t)+v(x,t)\\
u_t & = & w\!\!\!\!\!\nearrow_{\!\!\!t} +v_t=\alpha^2
(w\!\!\!\!\!\nearrow_{\!\!\!xx}
  +v_{xx})+s\!\!\!\!\!\!\nearrow (x,t)\\
\Rightarrow v_t & = & \alpha^2 v_{xx}.
\end{eqnarray}
$$(ref-heat-distributed-24)

$$
\begin{eqnarray}
\begin{array}{lclclclll}
A&=&u(0,t)&=&w(0,t)+v(0,t)&=&v(0,t)&\mbox{since }&w(0,t)=0\\
B&=&u(0,t)&=&w(L,t)+v(L,t)&=&v(L,t)&\mbox{since
}&w(L,t)=0.\end{array}\nonumber
\end{eqnarray}
$$(ref-heat-distributed-25)

$$
\begin{equation}
f(x)=u(x,0)=w(x,0)+v(x,0)\Rightarrow v(x,0)=f(x)-w(x,0)
\end{equation}
$$(ref-heat-distributed-26)
thus $v(x,t)$ satisfies:

$$
\begin{eqnarray}
\left.\begin{array}{ll}
&v_t=\alpha^2 v_{xx}\\ \mbox{BC: }&v(0,t)=A\quad v(L,t)=B\\
\mbox{IC: }&v(x,0)=f(x)-w(x,0)\end{array}\right\}
\label{eq_homog_dirichlet}
\end{eqnarray}
$$(eq_homog_dirichlet)

Now the boundary value problem {eq}`eq_homog_dirichlet` was solved
in lecture 19 of the notes.

Therefore

$$
\begin{eqnarray}
u(x,t)& = & \left(\frac{B-A}{L}\right) x+A+\sum\limits_{n=1}^\infty
\{\rm\ e\}^{-\alpha^2 {\left(\frac{n\pi}{L}\right)}^2 t}\left\{
b_n+\int\limits_0^{t} \{\rm\ e\}^{\alpha^2{\left(\frac{n\pi}{L}\right)}^2
\tau} \hat{s}_n (\tau )\, dx\right\} \sin\left(\frac{n\pi
x}{L}\right)
\end{eqnarray}
$$(ref-heat-distributed-28)

where

$$
\begin{eqnarray}
b_n=\frac{2}{L}\int\limits_0^L \left\{ f(x)-w(x,0)-\left[\left(
B-A\right)\frac{x}{L} +A\right]\right\}\sin\left(\frac{n\pi
x}{L}\right)\, dx.
\end{eqnarray}
$$(ref-heat-distributed-29)
````

````{prf:example} A bar with an external heat source dependent on $x$ and $t$: $S(x,t)=xt$
:label: example-heat-distributed-3 

$$
\begin{eqnarray*}
u_t&=&\alpha^2 u_{xx}+xt\\
u(0,t)&=&A\quad u(L,t)=B\\
u(x,t)&=&f(x)\\
\end{eqnarray*}
$$(ref-heat-distributed-30)

Let $q(x)=\left(\frac{B-A}{L}\right) x+A $ and
$u(x,t)=q(x)+v(x,t)$ then

$$
\begin{eqnarray*}
v_t&=&\alpha^2 v_{xx}+xt\\
v(0,t)&=&0=v(L,t)\\
v(x,0)&=&f(x)-q(x).
\end{eqnarray*}
$$(ref-heat-distributed-31)

Expanding the source $S(x,t)$ in terms of the eigenfunctions
of the problem with homogeneous boundary conditions

$$
\begin{eqnarray*}
S(x,t)=xt&=&\sum\limits_{n=1}^\infty \hat{S}_n(t)\sin (\lambda_nx)\\
\hat{S}_n(t)&=&\frac{2}{L}\int\limits_0^L xt\sin\left(\frac{n\pi x}{L}\right)\, dx\\
&=&\frac{2t}{L}\left\{ -x\frac{\cos\left(\frac{n\pi
x}{L}\right)}{\left(\frac{n\pi}{L}\right)}
  \Bigg|_{0\!\!\!\nearrow}^L +\left(\frac{L}{n\pi}\right)\int\limits_0^L\cos
  \left(\frac{n\pi x}{L}\right)\, dx\right\}\\
\hat{S}_n(t)&=&\frac{2t}{L}\left\{
(-1)^{n+1}\left(\frac{L^2}{n\pi}\right)\right\}
=\left(\frac{2L}{n\pi}\right) (-1)^{n+1}t
\end{eqnarray*}
$$(ref-heat-distributed-32)

$$
\begin{eqnarray*}
0=v_t-\alpha^2 v_{xx}-xt &=&\sum\limits_{n=1}^\infty\left\{
   \dot{\hat{v}}_n+\alpha^2\lambda_n^2\hat{v}_n-\hat{S}_n(t)\right\}
   \sin x \lambda_nx?\\
\dot{\hat{v}}_n+\alpha^2\lambda_n^2 \hat{v}_n&=&
  \left(\frac{2L}{n\pi}\right) (-1)^{n+1}t\\
\left( e^{\alpha^2\lambda_n^2 t}\hat{v}_n\right)
  &=&\left(\frac{2L}{n\pi}\right)
  (-1)^{n+1}
  \left[\frac{te^{\alpha^2\lambda_n^2t}}{\alpha^2\lambda_n^2}-
  \frac{e^{\alpha^2\lambda_n^2t}}{(\alpha^2\lambda_n^2)^2}\right]_0^t +c_n\\
&=&\left(\frac{2L}{n\pi}\right) (-1)^{n+1}\left[
  \frac{te^{\alpha^2\lambda_n^2t}}{\alpha^2\lambda_n^2}-
  \frac{(e^{\alpha^2\lambda_n^2t}-1)}{\alpha^4\lambda_n^4}\right] +c_n\\
\hat{v}_n(t)&=&\left(\frac{2L}{n\pi}\right) (-1)^n\left[
  \frac{t(\alpha^2\lambda_n^2) +e^{-\alpha^2\lambda_n^2t}-1}
   {(\alpha^2\lambda_n^2)^2}\right] +c_ne^{-\alpha^2\lambda_n^2t}\\
v(x,t)&=&\sum\limits_{n=1}^\infty\left\{\left(\frac{2L}{n\pi}\right)
(-1)^n
  \left[\frac{t(\alpha^2\lambda_n^2)-1+e^{-\alpha^2\lambda_n^2t}}
   {(\alpha^4\lambda_n^4)}\right]
   +c_ne^{-\alpha^2\lambda_n^2t}\right\}\sin (\lambda_nx)\\
f(x)-q(x)&=&v(x,0)=\sum\limits_{n=1}^\infty \left(
  \frac{2L}{n\pi}\right) (-1)^n 0+c_n\sin (\lambda_nx)\\
c_n&=&\frac{2}{L}\int\limits_0^L\left[
f(x)-\left\{\left(\frac{B-A}{L}\right) x+A\right\}\right]\sin
(\lambda_nx)\, dx
\end{eqnarray*}
$$(ref-heat-distributed-33)
````
