
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

# Full Range Fourier Series

In this lecture we consider the Full Range Fourier Series for a
given function defined on an interval $[-L,L]$. Outside this
interval we see that the Fourier Series represents the periodic
extension of the function $f(x)$.

```{admonition} Key Concepts
Full Range Fourier Series; Periodic
Extension; Complex Fourier Series.
```

## Fourier Series

We consider the expansion of the function $f(x)$ of the form

$$
\begin{equation}
f(x)\sim\frac{a_0}{2}+\sum\limits_{n=1}^\infty a_n
\cos\left(\frac{n\pi x}{L}\right) +b_n\sin\left(\frac{n\pi
x}{L}\right) = S(x) \label{eq:FS}
\end{equation}
$$(ref-fourier-full-0)

where

$$
\begin{eqnarray}
a_n & = & \frac{1}{L}\int\limits_{-L}^Lf(x)\cos\left(\frac{n\pi
x}{L}\right)\, dx\qquad\frac{a_0}{2}=\frac{1}{2L}\int\limits_{-L}^L
f(x)\, dx =\mbox{ average value of $f$.}\nonumber\\
b_n & = & \frac{1}{L}\int\limits_{-L}^L f(x)\sin\left(\frac{n\pi
x}{L}\right)\, dx
\end{eqnarray}
$$(ref-fourier-full-1)

````{prf:observation}
:label: observation-fourier-full-0
1. Note that $\displaystyle\cos\left(\frac{n\pi}{L}(x+T)\right)
   =\cos\left(\frac{n\pi x}{L}\right)$ provided $
   \displaystyle\frac{n\pi T}{L}=2\pi, T =\displaystyle\frac{2L}{n}
   $ and similarly $\displaystyle\sin\left(\frac{n\pi}{L}(x+2L)\right) 
   =\sin\left(\frac{n\pi x}{L}\right)$.

   Thus each of the terms of the Fourier Series $S(x)$
   on the RHS of (\ref{eq:FS}) is a periodic function having a maximal
   period $2L$. As a result the function $S(x)$ is also periodic.

   How does this relate to $f(x)$ which may not be
   periodic?

2. The function $S(x)$ represented by the series is known as the
   periodic extension of $f$ on $[-L,L]$.

   If $f$ (or its periodic extension) is discontinuous at a point
   $x_0$ then $S(x)$ converges to the average value of $f$ across the
   discontinuity.

   $$
   \begin{eqnarray}
   S(x_0)=\frac{1}{2}\left\{ f(x_0^+)+f(x_0^-)\right\}
   \end{eqnarray}
   $$(ref-fourier-full-2)
````

````{prf:example}
:label: example-fourier-full-0
$$
\begin{eqnarray}f(x)=\left\{\begin{array}{lll} 0\qquad &-\pi <x<0\quad &L=\pi\\
x&\phantom{-}0\leq x\leq\pi &\end{array}\right.\end{eqnarray}
$$(ref-fourier-full-3)

$$
\begin{eqnarray}
a_0 & = & \frac{1}{\pi}\int\limits_{-\pi}^\pi f(x)\,
dx=\frac{1}{\pi}\int\limits_0^\pi x\, dx=\frac{\pi}{2}\\
a_n & = & \frac{1}{\pi}\int\limits_{-\pi}^\pi f(x)\cos (nx)\, dx \nonumber\\
    & = & \frac{1}{\pi}\int\limits_0^\pi x\cos (nx)\, dx \nonumber\\
    & = & \frac{1}{\pi}\left\{\left. x\frac{\sin
    (nx)}{n}\right|_0^\pi -\frac{1}{n}\int\limits_0^\pi 1.\sin
    (nx)\, dx\right\} \nonumber \\
    & = & \frac{1}{\pi}\left\{\frac{\pi\sin}{n}(n\pi
    )\!\!\!\!\!\!\!\!\!\!\!\nearrow +\left.\frac{1}{n^2}\cos (nx)\right|_0^\pi\right\} \nonumber\\
    & = & \frac{1}{\pi n^2}\big[ (-1)^n-1\big]
    \quad\begin{array}{crcrc}
    n&1&2&3&4\\ (-1)^n-1&-2&0&-2&0\end{array}\\
a_{2m+1} & = & -\frac{2}{\pi (2m+1)^2}\quad m=0,1,2,\ldots\\
b_n & = & \frac{1}{\pi}\int\limits_{-\pi}^\pi f(x)\sin (nx)\, dx \nonumber\\
    & = & \frac{1}{\pi}\int\limits_0^\pi x\sin (nx)\, dx \nonumber\\
    & = & \frac{1}{\pi}\left\{\left. -x\frac{\cos
    (nx)}{n}\right|_0^\pi +\frac{1}{n}\int\limits_0^\pi 1.\cos
    (nx)\, dx\right\} \nonumber \\
    & = & \frac{1}{\pi}\left\{ -\pi\frac{\cos (n\pi )}{n}+
    \frac{0.\cos 0}{n}+\left.\frac{1}{n^2}\sin
    (nx)\right|_0^\pi \right\} \nonumber\\
    & = & (-1)^{n+1}/n\\
f(x) & = & \frac{a_0}{2}+\sum\limits_{n=1}^\infty a_n\cos
    (nx)+b_n\sin (nx) \nonumber\\
    & = & \frac{\pi}{4}-\frac{2}{\pi}\sum\limits_{m=0}^\infty
    \frac{\cos\big[ (2m+1)x\big]}{(2m+1)^2}+\sum\limits_{n=1}^\infty
    (-1)^{n+1}\frac{\sin (nx)}{n}
\end{eqnarray}
$$(ref-fourier-full-4)

```{figure} ../img/fourier/x_10terms.png
:name: x_10terms
:align: center

Truncated Fourier Series approximation to $f(x)$ using 10 terms. Notice the periodic extension of the function that was sampled on $[-\pi,\pi]$
and the oscillations in the Fourier Series near the points of discontinuity. Also note that at the point of discontinuity $x=\pi$, $ S(\pi)=\frac{1}{2}\left\{ f(\pi^+)+f(\pi^-)\right\}$
```
````

## It can be useful to shift the interval of integration from $[-L,L]$ to $[c,c+2L]$

Since the periodic extension $f_\{\rm\ e\} (x)$ is periodic with period
$2L$ (as are the basis functions $\displaystyle\cos\left(\frac{n\pi x}{L}\right)$ and $\displaystyle\sin\left(\frac{n\pi x}{L}\right)$).

$$
\begin{eqnarray}
a_n & = &\frac{1}{L}\int\limits_{-L}^L f(x)\cos\left(\frac{n\pi
x}{L}\right)\, dx
   =\frac{1}{L}\int\limits_c^{c+2L}f_\{\rm\ e\} (x)\cos\left(\frac{n\pi x}{L}\right)\, dx\\
b_n & = &\frac{1}{L}\int\limits_{-L}^Lf(x)\sin\left(\frac{n\pi
x}{L}\right)\, dx=\frac{1}{L}\int\limits_c^{c+2L}f_\{\rm\ e\}
(x)\sin\left(\frac{n\pi x}{L}\right)\, dx.
\end{eqnarray}
$$(ref-fourier-full-5)

````{prf:example}
:label: example-fourier-full-1

From the previous example,

$$
\begin{eqnarray}f(x)=\left\{\begin{array}{ll}
0 &-\pi < x< 0\\
x &\phantom{-}0\leq x\leq\pi\end{array}\right.
\end{eqnarray}
$$(ref-fourier-full-6)

On $[\pi ,3\pi ]$

$$
\begin{eqnarray}f_\{\rm\ e\} (x)=\left\{\begin{array}{ll}
0 &\phantom{2}\pi <x<2\pi\\
x-2\pi &2\pi\leq x\leq 3\pi\end{array}\right.
\end{eqnarray}
$$(ref-fourier-full-7)

$$
\begin{eqnarray}\begin{array}{lcl}
a_n &= &\frac{1}{\pi}\int\limits_\pi ^{3\pi}f_\{\rm\ e\} (x)\cos (nx)\, dx\\
 &= &\frac{1}{\pi}\int\limits_{2\pi}^{3\pi} (x-2\pi )\cos (nx)\, dx\\
 &= &\frac{1}{\pi}\int\limits_0^\pi t\cos (nt)\, dt.\end{array}\hspace{.5in}
\begin{array}{ll}
t=x-2\pi &dx=dt\\
x=t+2\pi &x=\pi\Rightarrow t=-\pi\\
  &x=3\pi\Rightarrow t=\pi\\
\mbox{since }&\cos n(t+2\pi )=\cos n t\end{array}
\end{eqnarray}
$$(ref-fourier-full-8)
````

## Complex Form of Fourier Series

$$
\begin{eqnarray}
f(x) & = &
\frac{a_0}{2}+\sum\limits_{n=1}^\infty a_n\cos\left(\frac{n\pi x}{L}\right) +b_n\sin\left(\frac{n\pi x}{L}\right)\nonumber\\
\cos\left(\frac{n\pi x}{L}\right) & = &
\frac{\{\rm\ e\}^{i\left(\frac{n\pi x}{L}\right)}+\{\rm\ e\}^{-i\left(\frac{n\pi
x}{L}\right)}}{2};\quad
\sin\left(\frac{n\pi x}{L}\right) =\frac{\{\rm\ e\}^{i\left(\frac{n\pi x}{L}\right)}-\{\rm\ e\}^{-i\left(\frac{n\pi x}{L}\right)}}{2i}\nonumber\\
\mbox{Therefore }f(x) & = & \frac{a_0}{2}+\sum\limits_{n=1}^\infty
\frac{a_n}{2} \left\{\{\rm\ e\}^{i\left(\frac{n\pi
x}{L}\right)}+\{\rm\ e\}^{-i\left(\frac{n\pi x}{L}\right)}\right\}
   +\frac{b_n}{2i}\left\{ \{\rm\ e\}^{i\left(\frac{n\pi x}{L}\right)}-\{\rm\ e\}^{-i\left(\frac{n\pi
    x}{L}\right)}\right\}\nonumber\\
& = & \frac{a_0}{2}+\sum\limits_{n=1}^\infty
\left(\frac{a_n-ib_n}{2}\right)
    \{\rm\ e\}^{i\left(\frac{n\pi x}{L}\right)}+\left(\frac{a_n+ib_n}{2}\right)
    \{\rm\ e\}^{-i\left(\frac{n\pi x}{L}\right)}\\
& &\, \uparrow\hspace{.85in}\uparrow\hspace{1.5in}\uparrow\nonumber\\
& & \, c_0\hspace{.85in}c_n\hspace{1.2in}c_{-n}\nonumber\\
& = & \sum\limits_{n=-\infty}^\infty c_n\{\rm\ e\}^{i\left(\frac{n\pi
x}{L}\right)}\nonumber
\end{eqnarray}
$$(ref-fourier-full-9)

$$
\begin{eqnarray}
c_n =\frac{a_n-ib_n}{2} & = &\frac{1}{2L}\int\limits_{-L}^L f(x)
   \left\{\cos\left(\frac{n\pi x}{L}\right) - i\sin\left(\frac{n\pi x}{L}\right)\right\}\, dx\\
& = & \frac{1}{2L}\int\limits_{-L}^L f(x)\{\rm\ e\}^{-i\left(\frac{n\pi
x}{L}\right)}\, dx\qquad b_{-n}=-b_n
\end{eqnarray}
$$(ref-fourier-full-10)

Therefore

$$
\begin{eqnarray}
f(x)& = &\sum\limits_{n=-\infty }^\infty c_n\{\rm\ e\}^{i\left(\frac{n\pi x}{L}\right)}\\
c_n & = &\frac{1}{2L}\int\limits_{-L}^L f(x)\{\rm\ e\}^{-i\left(\frac{n\pi
x}{L}\right)}\, dx.
\end{eqnarray}
$$(ref-fourier-full-11)

````{prf:example}
:label: example-fourier-full-2
$$
\begin{eqnarray}
f(x) & = &\left\{\begin{array}{rl}-1&-\pi\leq x<0\\
1&\phantom{-}0<x<\pi
\end{array}
     \quad L=\pi\right.\\
c_n & = &\frac{1}{2\pi}\left\{ -\int\limits_{-L}^0\{\rm\ e\}^{-inx}\, dx
   +\int\limits_0^\pi \{\rm\ e\}^{-inx}\, dx\right\}\\
& = & \frac{1}{2\pi}\left\{ -
   \frac{\left. \{\rm\ e\}^{-inx}\right|_{-\pi}^0}{(-in)} +
   \frac{\left. \{\rm\ e\}^{-inx}\right|_0^\pi}{(-in)}\right\}\\
& = & \frac{i}{2\pi n}\left\{ -2+\{\rm\ e\}^{+in\pi}+\{\rm\ e\}^{-in\pi}\right\}
=\left\{\begin{array}{cl}0&n\mbox{ even}\\(2/i\pi n)&n\mbox{
odd}\end{array}\right.
\end{eqnarray}
$$(ref-fourier-full-12)

Therefore 

$$
\begin{eqnarray}f(x) & = & \sum\limits_{n=-\infty}^\infty
   \frac{2}{\pi i(2n+1)} \{\rm\ e\}^{i\big( (2n+1)x\big)}.
\end{eqnarray}
$$(ref-fourier-full-13)
