
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

# Circular Domains

In this lecture we consider the solution of Laplace's equations on
domains that have a boundary that has at least one boundary segment
that comprises a circular arc.

```{admonition} Key Concepts
Laplace's equation; Circular domains;
Pizza Slice-shaped regions, Dirichlet and Mixed BC.
```

Reference Section: Boyce and Di Prima Section 10.8

## General Analysis of Laplace's Equation on Circular Domains

### Laplacian in Polar Coordinates

For domains whose boundary comprises part of a circle, it is
convenient to transform to polar coordinates. For this purpose the
Laplacian is transformed from cartesian coordinates $(x,y)$ to polar
coordinates $(r,\theta)$ as follows:

$$
\begin{eqnarray*}
r^2&=&x^2+y^2\\
\theta &=&\tan^{-1}\left(\frac{y}{x}\right)
\end{eqnarray*}
$$(ref-laplace-circular-0)

Differentiating with respect to $x$ and $y$ we obtain:

$$
\begin{eqnarray*}
2rr_x&=&2x\quad r_x=\frac{x}{r}\quad r_y=\frac{y}{r}\\
\theta_x&=&\frac{(-{\frac{y}{x^2})}} {\big( 1+(\frac{y}{x})^2\big)}
=-\frac{y}{x^2+y^2}=-\frac{y}{r^2}\quad
\theta_y=\frac{x}{r^2}\\
u(x,y)&=&U(r,\theta )\\
u_x&=&U_rr_x+U_\theta \theta_x=U_r\big(\frac{x}{r}\big)
+U_\theta\left( -\frac{y}{r^2}\right)\\
u_y&=&U_rr_y+U_\theta \theta_y=U_r\left(\frac{y}{r}\right)
+U_\theta\left(\frac{x}{r^2}\right)\\
u_{xx}&=&(U_r)_x r_x+U_rr_{xx}+(U_\theta
)_x\theta_x+U_\theta\theta_{xx}\\
&=&U_{rr}r_x^2+U_{r\theta}\theta_xr_x+U_rr_{xx}+U_{r\theta}r_x
  \theta_x+U_{\theta\theta}\theta_x^2+U_\theta\theta_{xx}\\
r_{xx}&=&\frac{r-\left(\frac{x^2}{r}\right)}{r^2}=\frac{y^2}{r^3}\quad
\theta_{xx}=\frac{2yx}{r^4}\\
u_{xx}&=&U_{rr}\left(\frac{x^2}{r^2}\right) +2U_{r\theta}\left(
-\frac{y}{r^2}\right)\left(\frac{x}{r}\right)
+U_{\theta\theta}\frac{y^2}{r^4}+U_r\left(\frac{y^2}{r^3}\right)
+U_\theta\left(\frac{2xy}{r^4}\right)\\
u_{yy}&=&U_{rr}r_y^2+U_{r\theta}r_y\theta_y+U_rr_{yy}+U_{\theta
r}r_y\theta_y+U_{\theta\theta}\theta_y^2+U_\theta\theta_{yy}\\
&=&U_{rr}\left(\frac{y^2}{r^2}\right)
+2U_{r\theta}\left(\frac{x}{r^2}\right)\left(\frac{y}{r}\right)
+U_{\theta\theta}\left(\frac{x^2}{r^4}\right)
+U_r\left(\frac{x^2}{r^3}\right) +U_\theta\left(
-\frac{2xy}{r^4}\right)
\end{eqnarray*}
$$(ref-laplace-circular-1)

$$
\begin{eqnarray*}
\framebox{$\displaystyle
u_{xx}+u_{yy}=U_{rr}+\frac{1}{r}U_r+\frac{1}{r^2} U_{\theta\theta}$}
\end{eqnarray*}
$$(ref-laplace-circular-2)

### Introductory Remarks About Circular Domains

Recall the Laplacian in polar coordinates:

$$
\begin{eqnarray}
0=\Delta u=u_{xx}+u_{yy}=u_{rr}+\frac{1}{r}u_r+\frac{1}{r^2}u_{\theta\theta}\quad\begin{array}{l}r=(x^2+y^2)^{1/2}\\
\theta =\tan^{-1} (y/x)\end{array}.
\end{eqnarray}
$$(ref-laplace-circular-3)

Let

$$
\begin{eqnarray}
u(r,\theta )&=&R(r)\Theta (\theta )\\
\frac{r^2R^{\prime\prime}+rR^\prime}{R(r)} &=&
-\frac{\Theta^{\prime\prime}}{\Theta (\theta )}=\lambda^2
\end{eqnarray}
$$(ref-laplace-circular-4)

which leads to $r^2R^{\prime\prime}+rR^\prime -\lambda^2 R=0$ and
$\Theta^{\prime\prime}+\lambda^2\Theta =0$.

The $R$ Equation: $r^2R^{\prime\prime}+rR^\prime-\lambda^2 R=0$:

$\mathbf{\lambda =0}$: $r^2R^{\prime\prime}+rR^\prime =0$,
$R=r^\gamma\Rightarrow\gamma(\gamma -1)+\gamma =\gamma^2 =0\Rightarrow R(r)=C+D\ln r$

$\mathbf{\lambda\ne 0}$: $r^2R^{\prime\prime}+rR^\prime -\lambda^2
=0$, $R=r^\gamma\Rightarrow\gamma (\gamma -1)+\gamma -\lambda^2
=\gamma^2 -\lambda^2=0\Rightarrow R(r)=Cr^\lambda +Dr^{-\lambda}$.

The $\Theta$ Equation $\Theta^{\prime\prime}+\lambda^2\Theta =0$:
$\Theta^{\prime\prime}+\lambda^2 \Theta =0$,
$\Theta =A\cos\lambda\theta +B\sin\lambda\theta$,
$\Theta^\prime =-A\lambda\sin\lambda\theta +B\lambda\cos\lambda\theta$  

Different Boundary Conditions and Corresponding Eigenfunctions:

**(I)** $\Theta (0)=0=\Theta (\alpha )$ $\lambda_n=n\pi /\alpha , \ $
$n=1,2,\ldots, \ $ $\Theta_n(\theta )=\sin\lambda_n\theta$

**(II)** $\Theta^\prime (0)=0=\Theta^\prime (\alpha ), \ $ $\lambda_n=n\pi
/\alpha $ $n=0,1,2,\ldots , \ $ $\Theta_n (\theta )\in\{
1,\cos\lambda_n\theta\}$

**(III)** $\Theta (0)=0=\Theta^\prime (\alpha ), \ $ $\lambda_n=(2n-1)\pi
/2\alpha$ $n=1,2,\ldots , \ $ $\Theta_n(\theta )=\sin\lambda_n \theta$

**(IV)** $\Theta^\prime (0)=0=\Theta (\alpha ), \ $ $\lambda_n=(2n-1)\pi
/2\alpha , \ $ $n=1,2,\ldots, \ $ $\Theta_n(\theta )=\cos\lambda_n\theta$

**(V)** $\displaystyle\left.\begin{array}{lcl}\Theta (-\pi )&=&\Theta (\pi
)\\ \Theta^\prime (-\pi )&=&\Theta^\prime (\pi )\end{array}\right\}$
$\lambda_n=n, \ $ $n=0,1,2,\ldots, \ $ $\Theta_n(\theta )\in\{
1,\cos\lambda_n\theta ,\sin\lambda_n\theta\}$.

The most general solution is thus of the form

$$
\begin{eqnarray}
u(r,\theta )=\left\{ A_0+\al_0 \ln r\right\}\cdot 1 &+&
\sum\limits_{n=1}^\infty \left\{ A_nr^{\lambda_n}+\al_n
r^{-\lambda_n}\right\}\cos\lambda_n\theta\\
&+&\sum\limits_{n=1}^\infty \left\{
B_nr^{\lambda_n}+\beta_nr^{-\lambda_n}\right\}\sin\lambda_n\theta .
\end{eqnarray}
$$(ref-laplace-circular-5)

````{prf:observation}
:label: observation-laplace-circular-0
- For problems that include the origin, the condition
$|u|<\infty$ as $r\rightarrow 0$ dictates that $\al_0=0$, $\al_n=0$ and
$\beta_n=0$.
- For problems that involve infinite domains the condition
$|u|<\infty$ as $r\rightarrow\infty$ dictates that $A_n=0$ and $B_n=0$.
- The values of $\lambda_n$ and the corresponding eigenfunctions
depend on the boundary conditions (I)--(V) that apply.
````

### Wedge Problems

````{prf:example}
:label: example-laplace-circular-0 Wedge with Homogeneous BC on $\theta =0$, $\theta =\alpha<2\pi$  
:label: EGHomogeneousWedge

```{figure} ../img/laplace/dirichlet_pizza.png
:name: dirichlet_pizza
:align: center

Homogeneous Dirichlet Boundary conditions on a wedge shaped
domain {eq}`eqLaplaceDirichletWedge`.
```

$$
\begin{eqnarray}
u_{rr}+\frac{1}{r}u_r+\frac{1}{r^2}u_{\theta\theta}& = & 0
   \qquad 0<r<a,\quad 0<\theta <\alpha\\
u(r,0)=0\quad u(r,\alpha ) & = & 0, \ u(r,\theta )\mbox{ bounded as
}r\rightarrow 0,\  u(a,\theta )=f(\theta ) \label{eqLaplaceDirichletWedge}
\end{eqnarray}
$$(eqLaplaceDirichletWedge)

Let $u(r,\theta )=R( r ) \cdot\Theta (\theta )$.

$$
\begin{eqnarray}
r^2\frac{(R^{\prime\prime}+\frac{1}{r}R^\prime )}{R}
   = -\frac{\Theta^{\prime\prime}(\theta )}{\Theta (\theta )}
   =\lambda^2\Rightarrow\begin{array}{rcl}
r^2R^{\prime\prime} +rR^\prime -\lambda^2 R&=&0\mbox{ Euler Eq.} \nonumber \\
\Theta^{\prime\prime}+\lambda^2\Theta&=&0\end{array} \nonumber
\end{eqnarray}
$$(ref-laplace-circular-7)

$u(r,0)=R( r )\Theta (0)=0\Rightarrow\Theta (0)=0$;
$u(r,\alpha )=R( r )\Theta (\alpha )=0\Rightarrow\Theta (\alpha )=0$

$$
\begin{eqnarray}
\begin{array}{l}\mbox{Eigenvalue}\\ \mbox{Problem}\end{array}
   \left\{\begin{array}{l}\Theta^{\prime\prime} +\lambda^2 \Theta =0\\
   \Theta (0)=0=\Theta (\alpha )\end{array}\quad\begin{array}{l}
   \Theta =A\cos\lambda \theta +B\sin (\lambda\theta )\\
   \Theta (0)=A=0\quad \Theta (\alpha )=B\sin (\lambda\alpha )=0\end{array}\right.
\end{eqnarray}
$$(ref-laplace-circular-8)

Therefore

$$
\begin{eqnarray}
\lambda_n =(n\pi /\alpha )\quad n=1,2, \ldots\quad\Theta_n =
   \sin\left(\frac{n\pi\theta}{\alpha}\right) .
\end{eqnarray}
$$(ref-laplace-circular-9)

To solve the Euler Eq. let $R( r )=r^\gamma$, $R^\prime =\gamma
r^{\gamma -1}$, $R^{\prime\prime}=\gamma (\gamma -1)r^{\gamma -2}$. Therefore

$$
\begin{eqnarray}
\gamma (\gamma -1)+\gamma -\lambda^2 =\gamma^2 -\lambda^2 =0
   \Rightarrow\gamma =\pm\lambda .
\end{eqnarray}
$$(ref-laplace-circular-10)

Therefore

$$
\begin{eqnarray}
R( r )=c_1r^\lambda +c_2r^{-\lambda}.
\end{eqnarray}
$$(ref-laplace-circular-11)

Now since $u(r,\theta )<\infty$ as $r\rightarrow 0$ we require $c_2=0$. Therefore

$$
\begin{eqnarray}
u(r,\theta ) & = & \sum\limits_{n=1}^\infty c_n
   r^{\left(\frac{n\pi}{\alpha}\right)} \sin
   \left(\frac{n\pi\theta}{\alpha}\right)\\
u(a,\theta )=f(\theta ) & = & \sum\limits_{n=1}^\infty
   \left\{ c_na^{\left(\frac{n\pi}{\alpha}\right)}\right\}
   \sin\left(\frac{n\pi\theta}{\alpha}\right) .
\end{eqnarray}
$$(ref-laplace-circular-12)

This is just a Fourier Sine Series for $f(\theta )$: Therefore

$$
\begin{eqnarray}
c_n a^{\left(\frac{n\pi}{\alpha}\right)} & = & \frac{2}{\alpha}
   \int\limits_0^\alpha f(\theta )\sin\left(\frac{n\pi\theta}{\alpha}\right)\, d\theta\\
c_n & = & \frac{2}{\alpha} a^{-\left(\frac{n\pi}{\alpha}\right)}
   \int\limits_0^\alpha f(\theta )\sin\left(\frac{n\pi\theta }{\alpha}\right)\, d\theta .
\end{eqnarray}
$$(ref-laplace-circular-13)

Therefore

$$
\begin{eqnarray}
u(x,\theta ) = \sum\limits_{n=1}^\infty c_n
r^{\left(\frac{n\pi}{\alpha}\right)}
   \sin\left(\frac{n\pi\theta}{\alpha}\right) .
\end{eqnarray}
$$(ref-laplace-circular-14)
````

````{prf:example}
:label: example-laplace-circular-1 A Wedge with Inhomogeneous BC
```{figure} ../img/laplace/inhomo_pizza.png
:name: inhomo_pizza
:align: center

Inhomogeneous Dirichlet Boundary conditions on a wedge
shaped domain {eq}`eqLaplaceInhomDirichletWedge`.
```

$$
\begin{eqnarray}
u_{rr}+\frac{1}{r}u_r+\frac{1}{r^2}u_{\theta\theta}
   =0\qquad 0<r<a,\quad 0<\theta <\alpha
\end{eqnarray}
$$(ref-laplace-circular-15)

$$
\begin{eqnarray}\begin{array}{l}
u(r,0)=u_0,\quad u(r,\alpha )=u_1,\quad u(r,\theta )<\infty\mbox{ as
}r\rightarrow 0,\quad u(a,\theta )=f(\theta )
\label{eqLaplaceInhomDirichletWedge}
\end{array}
\end{eqnarray}
$$(eqLaplaceInhomDirichletWedge)

Let us look for the simplest function of $\theta$ only that satisfies the
inhomogeneous BC of the from: $w(\theta
)=(u_1-u_0)\displaystyle\frac{\theta}{\alpha}+u_0$. Note that
$w_{\theta\theta}=0$ and that $w(0)=u_0$ and $w(\alpha )=u_1$. Then let
$u(r,\theta )=w(\theta )+v(r,\theta )$.

$$
\begin{eqnarray}
\left.\begin{array}{c}
u_{rr}+\displaystyle\frac{1}{r}u_r+\displaystyle\frac{1}{r^2}u_{\theta\theta}
   = v_{rr}+\displaystyle\frac{1}{r}v_r+\displaystyle\frac{1}{r^2}v_{\theta\theta} =0\\
v(r,0)=0\quad v(r,\alpha )=0\\
v(a,\theta )=f(\theta )-w(\theta
)\end{array}\right\}\begin{array}{l}
\mbox{Essentially the problem}\\
\mbox{solved in Example {prf:ref}`EGHomogeneousWedge`\end{array}
\end{eqnarray}
$$(ref-laplace-circular-17)

The solution is

$$
\begin{eqnarray}
u(r,\theta
)=(u_1-u_0)\frac{\theta}{\alpha}+u_0+\sum\limits_{n=1}^\infty
c_nr^{\left(\frac{n\pi}{\alpha}\right)}\sin\left(\frac{n\pi\theta}{\alpha}\right)
\end{eqnarray}
$$(ref-laplace-circular-18)

where

$$
\begin{eqnarray}
c_n=\frac{2}{\alpha } a^{-\left(\frac{n\pi}{\alpha}\right)}
   \int\limits_0^\infty \big[ f(\theta )-w(\theta )\big]\sin
   \left(\frac{n\pi\theta}{\alpha}\right)\, d\theta .
\end{eqnarray}
$$(ref-laplace-circular-19)
````

````{prf:example}
:label: example-laplace-circular-2 A Wedge with Insulating BC on $\theta =0$ and $\theta =\alpha<2\pi$
```{figure} ../img/laplace/mixed_pizza.png
:name: mixed_pizza
:align: center

Mixed Boundary conditions on a wedge shaped domain {eq}`eqLaplaceMixedWedge`.
```

$$
\begin{eqnarray}\begin{array}{l}
u_{rr}+\displaystyle\frac{1}{r}u_r+\displaystyle\frac{1}{r^2}u_{\theta\theta}=0\\
u_\theta (r,0)=0, \quad u_\theta (r,\alpha )=0, \  u(a,\theta
)=f(\theta )\label{eqLaplaceMixedWedge}\end{array}
\end{eqnarray}
$$(eqLaplaceMixedWedge)

Let

$$
\begin{eqnarray}
u(r,\theta )=R( r )\Theta (\theta )\Rightarrow r^2\left(
R^{\prime\prime} +\frac{1}{r}R^\prime \right) /R( r
)=-\Theta^{\prime\prime}/\Theta =\lambda^2
\end{eqnarray}
$$(ref-laplace-circular-21)

$\mathbf{\Theta \mbox{ equation}\rangle}$

$$
\begin{eqnarray}
& &\left.\begin{array}{l}
\Theta^{\prime\prime}+\lambda^2\Theta =0\\
\Theta^\prime (0)=0=\Theta^\prime (\alpha )\end{array}\right\}
\begin{array}{l}
\Theta (\theta )=A\cos\lambda\theta +B\sin (\lambda\theta )\\
\Theta^\prime (0)=B\lambda =0\;\lambda =0\mbox{ or }B=0;\end{array}\\
& &\hspace{1in}\begin{array}{l}
\Theta^\prime (\theta )=-A\lambda\sin (\lambda\theta )+B\lambda\cos (\lambda\theta )\\
\Theta^\prime (\alpha )=-A\lambda\sin (\lambda\alpha )=0\;
\lambda_n=\frac{n\pi}{\alpha};\;
    n=0,1,\ldots\end{array}
\end{eqnarray}
$$(ref-laplace-circular-22)

$\mathbf{R \mbox{ equation} \rangle}$
$r^2R_n^{\prime\prime}
   +rR_n^\prime -\lambda_n^2R_n=0$.

$\mathbf{n=0}$: $rR_0^{\prime\prime}
  +R_0^\prime ={(rR_0^\prime )}^\prime =0\Rightarrow rR_0^\prime
   =d_0\Rightarrow R_0( r )=c_0+d_0\ln r$.

$\mathbf{n\geq 1}$: $r^2R_n^{\prime\prime}
   +rR_n^\prime -\lambda^2 R_n=0\Rightarrow R_n
   =c_nr^{\lambda_n}+D_nr^{-\lambda_n}$.

Since $u(r,\theta )<\infty$ (i.e. must be bounded) as $r\rightarrow 0$ we
require $d_0=0=D_n$. Therefore

$$
\begin{eqnarray}
u(r,\theta ) & = & \frac{c_0}{2}+\sum\limits_{n=1}^\infty c_nr^{\left(\frac{n\pi}{\alpha}\right)}\cos\left(\frac{n\pi\theta}{\alpha}\right)\\
f(\theta )=u(a,\theta ) & = & \frac{c_0}{2}+\sum\limits_{n=1}^\infty c_na^{\left(\frac{n\pi}{\alpha}\right)}\cos\left(\frac{n\pi\theta }{\alpha}\right)\\
c_0 & = & \frac{2}{\alpha}\int\limits_0^\alpha f(\theta )d\theta\quad
   c_n=\frac{2}{\alpha}a^{-\left(\frac{n\pi}{\alpha}\right)}\int\limits_0^\alpha
   f(\theta )\cos\left(\frac{n\pi\theta}{\alpha}\right)\, d\theta\phantom{\int\int\int}\\
u(r,\theta ) & = & \frac{c_0}{2}+\sum\limits_{n=1}^\infty
c_nr^{\left(\frac{n\pi}{\alpha}\right)}\cos\left(\frac{n\pi\theta}{\alpha}\right)
.
\end{eqnarray}
$$(ref-laplace-circular-23)
````

````{prf:example}
:label: example-laplace-circular-3 Mixed BC - a 'crack like' problem.
```{figure} ../img/laplace/semicircle.png
:name: semicircle
:align: center

Mixed crack-like boundary conditions on a circular
domain as prescribed in {eq}`eqMixedCrack`.
```

$$
\begin{eqnarray}
\Delta u=u_{rr}+\frac{1}{r}u_r+\frac{1}{r^2}u_{\theta\theta}=0
\end{eqnarray}
$$(ref-laplace-circular-24)

subject to

$$
\begin{eqnarray}
u(r,0) & = & 0\quad \frac{\partial u}{\partial\theta}(r,\pi )=0 \label{eqMixedCrack}\\
u(a,\theta ) & = & f(\theta ). \nonumber
\end{eqnarray}
$$(eqMixedCrack)

Let $u(r,\theta )=R( r )\Theta (\theta )$.

$$
\begin{eqnarray}
r^2\frac{\left( R^{\prime\prime}+\frac{1}{r}R^\prime\right)}{R}
   = -\frac{\Theta^{\prime\prime}(\theta )}{\Theta (\theta )}=\lambda^2
\end{eqnarray}
$$(ref-laplace-circular-26)

$\mathbf{\Theta \mbox{ equation}\rangle}$

$$
\begin{eqnarray}\begin{array}{l}
\Theta^{\prime\prime}+\lambda^2 \Theta =0\\
\\
\Theta (0)=0\quad \Theta^\prime (\pi )=0\end{array}
\begin{array}{l}
\Theta =A\cos\lambda\theta +B\sin\lambda\theta\quad \Theta^\prime
   = - A\lambda\sin\lambda\theta +B\lambda\cos\lambda\theta\\
\Theta (0)=A=0\quad\Theta^\prime (\pi )=B\lambda\cos (\lambda\pi )
   =0\Rightarrow\pi\la_1 =\displaystyle\frac{\pi}{2}, \displaystyle\frac{3\pi}{2},
   \ldots\end{array}
\end{eqnarray}
$$(ref-laplace-circular-27)

or $\lambda_n =(2n+1)\displaystyle\frac{1}{2}$ $n=0,1,\ldots $ $\lambda\ne 0$ as this
would be trivial.

$\mathbf{R \mbox{ equation}\rangle}$
$r^2R^{\prime\prime}+rR^\prime -\lambda^2 R=0$ $R( r
)=r^\gamma\Rightarrow\gamma^2 -\lambda^2 =0$ $\gamma =\pm\lambda$. Therefore

$$
\begin{eqnarray}
u_n(r,\theta )=\left( c_nr^{\lambda_n}+d_nr^{-\lambda_n}\right)\sin\lambda_n
\theta .
\end{eqnarray}
$$(ref-laplace-circular-28)

Since $u$ should be bounded as $r\rightarrow 0$ we conclude that $d_n=0$.
The general solution is thus

$$
\begin{eqnarray}
u(r,\theta ) & = & \sum\limits_{n=0}^\infty c_nr^{(2n+1)/2}\sin\left(\frac{(2n+1)}{2}\theta\right)\\
f(\theta ) & = & u(a,\theta )=\sum\limits_{n=0}^\infty
c_na^{(2n+1)/2}\sin\left(\left(\frac{2n+1}{2}\right)\theta\right) .
\end{eqnarray}
$$(ref-laplace-circular-29)

Check orthogonality

$$
\begin{eqnarray}
\int\limits_0^\pi\sin\left(\left(\frac{2m+1}{2}\right)\theta\right)\sin\left(\left(\frac{2n+1}{2}\right)\theta\right)\,
d\theta =\left\{\begin{array}{ll}0&m\ne n\\ \pi
/2&m=n\end{array}\right. .
\end{eqnarray}
$$(ref-laplace-circular-30)

Therefore

$$
\begin{eqnarray}
c_n & = & \frac{2a^{-\left(
n+\frac{1}{2}\right)}}{\pi}\int\limits_0^\pi f(\theta
)\sin\left(\left( n+\frac{1}{2}\right)\theta\right)\, d\theta\\
u(r,\theta ) & = & \sum\limits_{n=0}^\infty c_nr^{\left(
  n+\frac{1}{2}\right)}\sin\left(\left( n+\frac{1}{2}\right)\theta\right)
\end{eqnarray}
$$(ref-laplace-circular-31)
````
