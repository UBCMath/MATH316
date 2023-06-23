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

# Heat Equation On A Circular Ring - Full Fourier Series

In this lecture we use separation of variables to solve the heat
equation subject on a thin circular ring with periodic boundary
conditions. In this case we reduce the problem to expanding the
initial condition function $f(x)$ in an infinite series of both
cosine and sine functions, which we refer to as the Full Range
Fourier Series.

```{admonition} Key Concepts
Heat Equation; Periodic Boundary
Conditions; separation of variables; Full Fourier Series.
```

```{figure} ../img/heat/circle.png
:name: circle
:alt: circle
:align: center

Consider a thin conducting ring with thermal conductivity
$\alpha^2$ that has a given initial temperature distribution
```

__Physical Interpretation:__ Consider a thin circular
wire in which there is no radial temperature dependence, i.e.,
$u(r,\theta)=u(\theta)$ so that $\displaystyle\frac{\partial u}{\partial
r}=0$. In this case the polar Laplacian reduces to

$$
\begin{eqnarray}
\Delta u & = & \frac{\partial^2u}{\partial
r^2}+\frac{1}{r}\frac{\partial u}{\partial
r}+\frac{1}{r^2}\frac{\partial^2 u}{\partial\theta^2}\nonumber\\
& = & \frac{\partial^2 u}{\partial (r\theta )^2}
\end{eqnarray}
$$(ref-heat-circular-0)

and if we let $x=r\theta$ then $\frac{\partial^2 u}{\partial
(r\theta )^2}=u_{xx}$. In this case the heat distribution in the
ring is determined by the following initial value problem with
periodic boundary conditions

$$
\begin{eqnarray}
u_t&=&\alpha^2u_{xx} \label{eq:heat}
\end{eqnarray}
$$(heat)

$$
\begin{eqnarray*}\begin{array}{llcl}
\mbox{BC:}&   & &\left.\begin{array}{lcl}u(-L,t)&=&u(L,t) \label{eq:ringheat}\\
                    \displaystyle\frac{\partial u}{\partial x}(-L,t)&=&\displaystyle\frac{\partial
                    u}{\partial
                    x}(L,t)\end{array}\right\}\quad\mbox{Periodic BC}\\
&&&\\
\mbox{IC:}&   & &u(x,0)=f(x)
\end{array}
\end{eqnarray*}
$$(ringheat)

Assume $\displaystyle u(x,t)=X(x)T(t)$.

As before:

$$
\displaystyle\frac{X^{\prime\prime}(x)}{X(x)}=\displaystyle\frac{\dot{T}(t)}{\alpha^2
T(t)}=-\lambda^2
$$(ref-heat-circular-3)

IVP:

$$
\displaystyle\frac{\dot{T}(t)}{\alpha^2 T(t)}=-\lambda^2\Rightarrow
T(t)=c\{\rm\ e\}^{-\lambda^2 t}
$$(ref-heat-circular-4)

$$
\displaystyle r=\frac{2L}{2\pi}=\frac{L}{\pi}=\mbox{ Constant}
$$(ref-heat-circular-5)

BVP:

$$
\left.\begin{array}{lcl} X^{\prime\prime}+\lambda^2X&=&0\\
X(-L)&=&X(L)\\ X^\prime (-L)&=&X^\prime (L)\end{array}\right\}
\quad\begin{array}{l}\mbox{Eigenvalue Problem}\\
\mbox{look for $\lambda$ such that}\\
\mbox{nontrivial $x$ can be found.}\end{array}
$$(ref-heat-circular-6)

$$
\begin{eqnarray}
X(x)& = & A\cos (\lambda x)+B\sin (\lambda x)\nonumber\\
X(-L)& = & A\cos (\lambda L)-B\sin (\lambda L)=A\cos (\lambda L)+B\sin (\lambda
L)=X(L)\nonumber\\
& &\quad\mbox{therefore }2B\sin (\lambda L)=0\nonumber\\
X^\prime (x) & = & -A\lambda \sin\lambda x+B\lambda\cos (\lambda x)\\
X^\prime (-L)& = & +A\lambda\sin (\lambda L)+B\lambda\cos (\lambda L)=-A\lambda\sin (\lambda
L)+B\lambda\cos (\lambda L)=X^\prime (L)\nonumber\\
& &\quad\mbox{therefore }2A\lambda\sin (\lambda L)=0\nonumber\\
\mbox{Therefore }\la_nL & = &(n\pi )\quad n=0,1,\ldots .\nonumber
\end{eqnarray}
$$(ref-heat-circular-7)

Solutions to {eq}`heat` that satisfy the BC are thus of the form

$$
\begin{equation}
u_n(x,t)=\{\rm\ e\}^{-{\left(\frac{n\pi}{L}\right)}^2\alpha^2t}
  \left\{A_n\cos\left(\frac{n\pi x}{L}\right) +B_n\sin
  \left(\frac{n\pi x}{L}\right)\right\} .
\end{equation}
$$(ref-heat-circular-8)

Superposition of all these solutions yields the general solution

$$
\begin{equation}
u(x,t)=A_0+\sum\limits_{n=1}^\infty\left\{ A_n\cos\left(\frac{n\pi
x}{L}\right) +B_n\sin\left(\frac{n\pi
x}{L}\right)\right\}\{\rm\ e\}^{-{\left(\frac{n\pi}{L}\right)}^2\alpha^2t}.
\end{equation}
$$(ref-heat-circular-9)

In order to match the IC we have

$$
\begin{equation}
f(x)=u(x,0)=A_0+\sum\limits_{n=1}^\infty A_n\cos\left(\frac{n\pi
x}{L}\right) +B_n\sin\left(\frac{n\pi x}{L}\right)
\label{eq:FourierSeries}.
\end{equation}
$$(FourierSeries)

As before we obtain expressions for the $A_n$ and $B_n$ by
projecting $f(x)$ onto $\displaystyle\sin\left(\frac{n\pi x}{L}\right)$ and
$\displaystyle\cos\left(\frac{n\pi x}{L}\right)$.

$$
\begin{eqnarray}
\int\limits_{-L}^Lf(x)\left\{\begin{array}{c} \sin\left(\frac{m\pi
x}{L}\right)\\
\cos\left(\frac{m\pi x}{L}\right)\end{array}\right\}\, dx & = &
A_0\int\limits_{-L}^L\left\{\begin{array}{c}\sin\left(\frac{m\pi
x}{L}\right)\\
\cos\left(\frac{m\pi x}{L}\right)\end{array}\right\}\, dx \label{eq:FourierProjection}\\
& &\quad +\, \sum\limits_{n=1}^\infty
A_n\int\limits_{-L}^L\cos\left(\frac{n\pi x}{L}\right)
\left\{\begin{array}{c}\sin\left(\frac{m\pi x}{L}\right)\\
\cos\left(\frac{m\pi x}{L}\right)\end{array}\right\}\, dx\nonumber\\
& &\quad +\, \sum\limits_{n=1}^\infty B_n\int\limits_{-L}^L
\sin\left(\frac{n\pi x}{L}\right)\left\{\begin{array}{c}
\sin\left(\frac{m\pi x}{L}\right)\\
\cos\left(\frac{m\pi x}{L}\right)\end{array}\right\}\, dx.\nonumber
\end{eqnarray}
$$(FourierProjection)

As in the previous example we use the orthogonality relations:

$$
\begin{eqnarray}
\int\limits_{-L}^L\sin\left(\frac{m\pi
x}{L}\right)\sin\left(\frac{n\pi x}{L}\right)\, dx & = &
  L\delta_{mn}\nonumber\\
\int\limits_{-L}^L\cos\left(\frac{m\pi x}{L}\right)\cos
  \left(\frac{n\pi x}{L}\right)\, dx & = &
  L\delta_{mn}\quad m\mbox{ and }n\not= 0\\
& = & 2L\quad m=n=0\nonumber\\
\int\limits_{-L}^L\sin\left(\frac{m\pi x}{L}\right)\cos
  \left(\frac{n\pi x}{L}\right)\, dx & = &
    0\quad\forall m,n.\nonumber
\end{eqnarray}
$$(ref-heat-circular-12)

Plugging these orthogonality conditions into {eq}`FourierProjection` we obtain

$$
\begin{eqnarray}\left.\begin{array}{lcl}
A_0 &= &\displaystyle\frac{1}{2L}\int\limits_{-L}^L f(x)\, dx=\mbox{ average
value of $f(x)$ on $[-L,L]$}\\
A_n &= &\displaystyle\frac{1}{L}\int\limits_{-L}^Lf(x)\cos\left(\frac{n\pi
x}{L}\right)\, dx\mbox{ and
}B_n=\frac{1}{L}\int\limits_{-L}^Lf(x)\sin\left(\frac{n\pi
x}{L}\right)\, dx.\end{array}\right\} \label{eq:FourierCoefficients}
\end{eqnarray}
$$(FourierCoefficients)

````{prf:observation}
:label: observation-heat-circular-0
1. {eq}`FourierSeries` and {eq}`FourierCoefficients`
 represent the full Fourier Series Expansion
for $f(x)$ on the interval $[-L,L]$.

2. By defining $\displaystyle a_n=\frac{1}{L}\int\limits_{-L}^L
f(x)\cos\left(\frac{n\pi x}{L}\right)\, dx=\left\{\begin{array}{c}
  2A_0\\ \\ A_n\end{array}\right.$ and $b_n=B_n$ the Fourier Series
{eq}`FourierSeries` is often written in the form

$$
\begin{eqnarray}
f(x)=\frac{a_0}{2}+\sum\limits_{n=1}^\infty a_n\cos\left(\frac{n\pi
x}{L}\right) +b_n\sin\left(\frac{n\pi x}{L}\right) .
\end{eqnarray}
$$(ref-heat-circular-14)
````
