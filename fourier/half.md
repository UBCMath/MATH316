
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

# Half Range Fourier Series: Even And Odd Functions

In this lecture we consider the Fourier Expansions for Even and Odd
functions, which give rise to cosine and sine half range Fourier
Expansions. If we are only given values of a function $f(x)$ over
half of the range $[0,L]$, we can define two different extensions of
$f$ to the full range $[-L,L]$, which yield distinct Fourier
Expansions. The even extension gives rise to a half range cosine
series, while the odd extension gives rise to a half range sine
series.

```{admonition} Key Concepts
Even and Odd Functions; Half Range
Fourier Expansions; Even and Odd Extensions
```

## Even and Odd Functions

__Even:__ $f(-x)=f(x)$

__Odd:__ $f(-x)=-f(x)$

### Integrals of Even and Odd Functions

$$
\begin{eqnarray}
\int\limits_{-L}^Lf(x)\, dx & = & \int\limits_{-L}^0f(x)\, dx+\int\limits_0^Lf(x)\, dx\\
& = & \int\limits_0^L\big[ f(-x)+f(x)\big]\, dx\\
& = & \left\{\begin{array}{ll}2\int\limits_0^Lf(x)\, dx &f\mbox{ even}\\
     0 &f\mbox{ odd}.\end{array} \right.
\end{eqnarray}
$$(ref-fourier-half-0)

```{note}
Let $E(x)$ represent an even function and $O(x)$ an
odd function.
- If $f(x)=E(x)\cdot O(x)$ then $f(-x)=E(-x)O(-x)=-E(x)O(x)=-f(x)\Rightarrow f$ is odd.
- $E_1(x)\cdot E_2(x)\rightarrow$ even.
- $O_1(x)\cdot O_2(x)\rightarrow$ even.
- Any function can be expressed as a sum of an even part and an odd part:

$$
\begin{eqnarray}
f(x)=\frac{1}{2}\underbrace{\big[ f(x)+f(-x)\big]}_{\mbox{even
part}}
   + \frac{1}{2}\underbrace{\big[ f(x)-f(-x)\big]}_{\mbox{odd part}}.
\end{eqnarray}
$$(ref-fourier-half-1)

__Check:__ Let $\displaystyle E(x)=\frac{1}{2}\big[ f(x)+f(-x)\big]$.
Then $\displaystyle E(-x)=\frac{1}{2}\big[ f(-x)+f(x)\big] =E(x)$ even.

Similarly let

$$
\begin{eqnarray}
O(x) & = & \frac{1}{2}\big[ f(x)-f(-x)\big]\\
O(-x) & = & \frac{1}{2}\big[ f(-x)-f(x)\big] =-O(x)\mbox{ odd}.
\end{eqnarray}
$$(ref-fourier-half-2)
```

## Consequences of the Even/Odd Property for Fourier Series

1. Let $f(x)$ be Even-Cosine Series:

   $$
   \begin{eqnarray}
   a_n & = & \frac{1}{L}\int\limits_{-L}^L
   \underbrace{f(x)\cos}_{\mbox{even}}\left(\frac{n\pi x}{L}\right)\,
      dx =\frac{2}{L}\int\limits_0^L f(x)\cos\left(\frac{n\pi x}{L}\right)\, dx\\
   b_n & = & \frac{1}{L}\int\limits_{-L}^L\underbrace{f(x)
      \sin\left(\frac{n\pi x}{L}\right)}_{\mbox{odd}}\, dx=0.
   \end{eqnarray}
   $$(ref-fourier-half-3)

   Therefore

   $$
   \begin{eqnarray}
   f(x)=\frac{a_0}{2}+\sum\limits_{n=1}^\infty a_n\cos\left(\frac{n\pi
   x}{L}\right) ;\quad a_n=\frac{2}{L}\int\limits_0^L
   f(x)\cos\left(\frac{n\pi x}{L}\right)\, dx.
   \end{eqnarray}
   $$(ref-fourier-half-4)

2. Let $f(x)$ be Odd-Sine Series:

   $$
   \begin{eqnarray}
   a_n & = & \frac{1}{L}\int\limits_{-L}^L\underbrace{f(x)\cos\left(\frac{n\pi x}{L}\right)}_{\mbox{odd}}\, dx=0\\
   b_n & = &
   \frac{1}{L}\int\limits_{-L}^L\underbrace{f(x)\sin\left(\frac{n\pi
   x}{L}\right)}_{\mbox{even}}\, dx
   =\frac{2}{L}\int\limits_0^Lf(x)\sin\left(\frac{n\pi x}{L}\right)\,
   dx \nonumber
   \end{eqnarray}
   $$(ref-fourier-half-5)

   Therefore

   $$
   \begin{eqnarray} f(x) & = &\sum\limits_{n=1}^\infty
   b_n\sin\left(\frac{n\pi x}{L}\right) ;\quad
   b_n=\frac{2}{L}\int\limits_{-0}^L f(x)\sin\left(\frac{n\pi
   x}{L}\right)\, dx. \nonumber
   \end{eqnarray}
   $$(ref-fourier-half-6)

3. Since any function can be written as the sum of an even and odd part,
   we can interpret the $\cos$ and $\sin$ series as even/odd:

   $$
   \begin{eqnarray}
   f(x) & = &
   \begin{array}{cc}
   \mbox{even} &\mbox{odd}\\
   \displaystyle\frac{1}{2}\big[ f(x)+f(-x)\big] &+\displaystyle\frac{1}{2}\big[ f(x)-f(-x)\big]\end{array}\\
   & = & \left\{ \frac{a_0}{2}+\sum\limits_{n=1}^\infty
   a_n\cos\left(\frac{n\pi x}{L}\right)\right\}
   +\left\{ \sum\limits_{n=1}^\infty b_n\sin\left(\frac{n\pi
   x}{L}\right)\right\} \nonumber
   \end{eqnarray}
   $$(ref-fourier-half-7)

   where

   $$
   \begin{eqnarray} a_n & = &
   \frac{2}{L}\int\limits_0^L\frac{1}{2}\big[ f(x)+f(-x)\big]\cos\left(
   \frac{n\pi x}{L}\right) \, dx
      =\frac{1}{L}\int\limits_{-L}^L f(x)\cos\left(\frac{n\pi x}{L}\right)\, dx \nonumber\\
   b_n & = & \frac{2}{L}\int\limits_0^L\frac{1}{2}\big[
   f(x)-f(-x)\big] \sin\left(\frac{n\pi x}{L}\right)\, dx
   =\frac{1}{L}\int\limits_{-L}^L f(x)\sin\left(\frac{n\pi x}{L}\right)
   \, dx. \nonumber
   \end{eqnarray}
   $$(ref-fourier-half-8)

## Half-Range Expansions

If we are given a function $f(x)$ on an interval $[0,L]$ and we want to
represent $f$ by a Fourier Series we have two choices - a Cosine Series
or a Sine Series.

__Cosine Series:__

$$
\begin{eqnarray}
f(x) & = & \frac{a_0}{2}+\sum\limits_{n=1}^\infty a_n\cos\left(\frac{n\pi x}{L}\right)\\
a_n & = & \frac{2}{L}\int\limits_0^L f(x)\cos\left(\frac{n\pi
x}{L}\right)\, dx.
\end{eqnarray}
$$(ref-fourier-half-9)

__Sine Series:__

$$
\begin{eqnarray}
f(x) & = & \sum\limits_{n=1}^\infty b_n\sin\left(\frac{n\pi x}{L}\right)\\
b_n & = & \frac{2}{L}\int\limits_0^L f(x)\sin\left(\frac{n\pi
x}{L}\right)\, dx.
\end{eqnarray}
$$(ref-fourier-half-10)

````{prf:example}
:label: example-fourier-half-0
Expand $f(x)=x$, $0<x<2$ in a half-range (a) Sine Series, (b) Cosine
Series.

**(a) Sine Series: (L=2)**

$$
\begin{eqnarray}
b_n & = & \frac{2}{L}\int\limits_0^L f(t)\sin\frac{n\pi}{\ell} t\, dt\\
& = & \int\limits_0^2 t\sin\frac{n\pi}{2} t\, dt\\
& = &
-\left.\frac{t\cos\frac{n\pi}{2}t}{\left(\frac{n\pi}{2}\right)}\right|_0^2
   +\frac{2}{n\pi}\int\limits_0^2\cos\frac{n\pi}{2}t\, dt\\
& = & -\frac{4}{n\pi}\cos (n\pi )+{\left(\frac{2}{n\pi}\right)}^2
   \left.\sin\left(\frac{n\pi}{2} t\right)\right|_0^2\\
& = & -\frac{4}{n\pi}(-1)^n
\end{eqnarray}
$$(ref-fourier-half-11)

Therefore

$$
\begin{eqnarray}
f(t) & = &\frac{4}{\pi}\sum\limits_{n=1}^\infty
\frac{(-1)^{n+1}}{n}\sin\left(\frac{n\pi}{2} t\right) .
\end{eqnarray}
$$(ref-fourier-half-12)

$$
\begin{eqnarray}
     f(1)&=&1=\frac{4}{\pi}\sum\limits_{n=1}^\infty
     \frac{(-1)^{n+1}}{n}\sin\left(\frac{n\pi}{2}\right)\\
\mbox{therefore }\
\frac{\pi}{4}&=&1-\frac{1}{3}+\frac{1}{5}-\frac{1}{7}+\cdots
\end{eqnarray}
$$(ref-fourier-half-13)

**(b) Cosine Series: (L=2)**

$$
\begin{eqnarray}
a_0 & = & \frac{2}{2}\int\limits_0^2 t\, dt=\left.\frac{t^2}{2}\right|_0^2 =2\\
a_n & = & \int\limits_0^2 t\cos\frac{n\pi}{2}t\,
dt=\left(\frac{2}{n\pi}\right) t\sin\!\!\!\!\!\!\!\!\nearrow
  \left.\frac{n\pi}{2} t\right|_0^2 -\left(\frac{2}{n\pi}\right) \int\limits_0^2\sin\frac{n\pi}{2} t\, dt \nonumber\\
& = & +\left. {\left(\frac{2}{n\pi}\right)}^2 \cos\frac{n\pi}{2}
t\right|_0^2
   =\frac{4}{n^2\pi^2}\left\{\cos n\pi -1\right\}
\end{eqnarray}
$$(ref-fourier-half-14)

Therefore

$$
\begin{eqnarray}
f(t) & = & 1+\frac{4}{\pi^2}\sum\limits_{n=1}^\infty \frac{\big[ (-1)^n-1\big]}{n^2} \cos\frac{n\pi}{2} t\\
& = & 1-\frac{8}{\pi^2}\sum\limits_{n=0}^\infty
\cos\frac{(2n+1)}{2}\pi t/(2n+1)^2.
\end{eqnarray}
$$(ref-fourier-half-15)

The cosine series converges faster than Sine Series.

$$
\begin{equation*}
f(2)=2=1+\frac{8}{\pi
^2}\sum\limits_{n=0}^\infty\frac{1}{(2n+1)^2},\quad \quad
\frac{\pi^2}{8}=1+\frac{1}{3^2}+\frac{1}{5^2}+\cdots
\end{equation*}
$$(ref-fourier-half-16)
````

````{prf:example}
:label: example-fourier-half-1 Periodic Extension
Assume that $f(x)=x$, $0<x<2$ represents one full period of the 
function so that $f(x+2)=f(x)$. $2L=2\Rightarrow
L=1$.

$$
\begin{eqnarray}
a_0 & = & \frac{1}{L}\int\limits_{-L}^L f(x)\, dx=\int\limits_{-1}^1 f(x)\, dx =\int\limits_0^2 x\, dx=\left.\frac{x^2}{2}\right|_0^2=2\\
&\phantom{=} & \hspace{1.5in}\mbox{since $f(x+2)=f(2)$.}
\end{eqnarray}
$$(ref-fourier-half-17)

$n\geq 1$:

$$
\begin{eqnarray}
a_n & = &\frac{1}{L}\int\limits_{-L}^L f(x)\cos\left(\frac{n\pi
x}{L}\right)\, dx
   = \int\limits_{-1}^1f(x)\cos (n\pi x)\, dx\quad L=1 \nonumber\\
& = & \int\limits_0^2 x\cos (n\pi x)\, dx \nonumber\\
& = & \left[\left.\frac{x\sin(n\pi
x)}{n\pi}\right|_{0\!\!\!\searrow}^{2\!\!\!\nearrow}
   -\left(\frac{1}{n\pi}\right)\int\limits_0^2\sin (n\pi x)\, dx\right]\nonumber \\
& = & \left.\frac{1}{(n\pi )^2}\cos (n\pi x)\right|_0^2
=\frac{1}{(n\pi )^2}
   \big[\cos (2n\pi )-1\big] =0 \\
b_n & = & \frac{1}{L}\int\limits_{-L}^L f(x)\sin\left(\frac{n\pi
x}{L}\right)\, dx
   =\int\limits_{-1}^1 f(x)\sin (n\pi x)\, dx \nonumber \\
& = & \int\limits_0^2 x \sin (n\pi x)\, dx =\left[
  -x\frac{\cos (n\pi x)}{n\pi}\Big|_{0\!\!\!\searrow}^{2}
  + \frac{1}{(n\pi )} \int\limits_0^2\cos (n\pi x)\, dx\right] \nonumber\\
& = & \frac{-2}{n\pi}+\left.\frac{\sin (n\pi x)}{(n\pi )^2}
  \right|_{0\!\!\!\searrow}^{2\!\!\!\nearrow} =\left(\frac{-2}{n\pi}\right)
\end{eqnarray}
$$(ref-fourier-half-18)

Therefore

$$
\begin{eqnarray}
f(x) & = & \frac{2}{2}-\frac{2}{\pi}\sum\limits_{n=1}^\infty \frac{\sin (n\pi x)}{n}\\
& = & 1  -\frac{2}{\pi}\sum\limits_{n=1}^\infty \frac{\sin (n\pi x)}{n} \nonumber\\
\end{eqnarray}
$$(ref-fourier-half-19)

```{figure} ../img/fourier/full_range.png
:name: full_range
:align: center

Left figure: Full Range Expansion $S_N(x)
= 1-\frac{2}{\pi}\sum\limits_{n=1}^{N=20} \frac{\sin (n\pi x)}{n}$.
Right figure: An odd function $S_N(x)-1
=-\frac{2}{\pi}\sum\limits_{n=1}^{N=20} \frac{\sin (n\pi x)}{n}$.
```
````
