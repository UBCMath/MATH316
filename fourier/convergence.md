
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

# Convergence of Fourier Series

In this lecture we state the fundamental convergence theorem for
Fourier Series, which assumes that the function $f(x)$ is piecewise
continuous. At points of discontinuity of $f(x)$ the Fourier
Approximation $S_N(x)$ takes on the average value
$\displaystyle\frac{1}{2}\big[ f(x+)+f(x-)\big]$ and exhibits the so-called
Gibbs Phenomenon in which the convergence is  pointwise but not
uniform. We explore the Gibbs phenomenon for a simple step function.

```{admonition} Key Concepts
Convergence of Fourier Series,
Piecewise continuous Functions, Gibbs Phenomenon.
```

- What conditions do we need to impose on $f$ to ensure that the Fourier Series converges to $f$.
- We consider piecewise continuous functions:

````{prf:theorem}
Let $f$ and $f^\prime$ be piecewise continuous functions on $[-L,L]$
and let $f$ be periodic with period $2L$, then $f$ has a Fourier
Series

$$
\begin{eqnarray}\begin{array}{llcl}
&f(x) &\sim &\frac{a_0}{2}+\sum\limits_{n=1}^\infty
a_n\cos\left(\frac{n\pi x}{L}\right)
   +b_n\sin\left(\frac{n\pi x}{L}\right) =S(x)\\
\mbox{where}&&&\\
&a_n &= &\frac{1}{L}\int\limits_{-L}^L f(x)\cos\left(\frac{n\pi
x}{L}\right)\, dx\mbox{ and }b_n=\frac{1}{L}\int\limits_{-L}^L
f(x)\sin\left(\frac{n\pi x}{L}\right)\, dx.\end{array}
\end{eqnarray}
$$(ref-fourier-convergence-0)

The Fourier Series converges to $f(x)$ at all points at which $f$ is
continuous and to $\displaystyle\frac{1}{2}\big[ f(x+)+f(x-)\big]$ at all
points at which $f$ is discontinuous.
````

- Thus a Fourier Series converges to the average value of the left and
  right limits at a point of discontinuity of the function $f(x)$.

## Illustration of the Gibbs Phenomenon - Nonuniform Convergence

- Near points of discontinuity truncated Fourier Series exhibit oscillations - overshoot.

```{figure} ../img/fourier/gibbs.png
:name: gibbs
:align: center

Fourier Series for a step function.
```

````{prf:example}
:label: example-fourier-convergence-0

Consider the half-range sine series expansion of

$$
\begin{equation}
f(x)=1\qquad\mbox{on }[0,\pi ].
\end{equation}
$$(ref-fourier-convergence-1)

$$
\begin{eqnarray}\begin{array}{rcll}
f(x)=1 &= &\sum\limits_{n=1}^\infty b_n\sin (nx) &\\
\mbox{where }b_n &= &\frac{2}{\pi}\int\limits_0^\pi\sin (nx)\, dx &=
\frac{2}{\pi}
   \left[ -\frac{\cos nx}{n}\right]_0^\pi =\frac{2}{\pi n}\big[ 1-(-1)^n\big]\\
&= & \left\{\begin{array}{lll}4/\pi n &n &\mbox{odd}\\
   0 &n&\mbox{even}\end{array}\right. &\\
\mbox{Therefore }f(x) &=
&\frac{4}{\pi}\sum\limits_{\stackrel{n=1}{n\mbox{ odd}}}
   ^\infty \frac{\sin (nx)}{n} &= \frac{4}{\pi }\sum\limits_{m=0}^\infty\frac{\sin (2m+1)x}{(2m+1)}.\end{array}
\end{eqnarray}
$$(ref-fourier-convergence-2)
````

```{note}
- $\displaystyle f(\pi /2)=1=\frac{4}{\pi}\sum\limits_{m=0}^\infty \frac{\sin\big[
  (2m+1)\pi /2\big]}{(2m+1)}=\frac{4}{\pi }\left\{ 1-\frac{1}{3}+\frac{1}
  {5}-\cdots\right\}$.
  Therefore $\displaystyle\frac{\pi}{4}=1-\frac{1}{3}+\frac{1}{5}-\cdots$.
- Recall the complex Fourier Series example for the function

$$
\begin{equation}
f(x)=\left\{
\begin{array}{c}
-1\text{ \ }-\pi \leq x<0 \\
1\text{ \ }\phantom{-}0<x<\pi
\end{array}
\right.   \label{function}
\end{equation}
$$(ref-fourier-convergence-3)

which turns out to be equivalent to the odd extension of
the above function represented by the half-range sine expansion,
which we can see from the following calculation

$$
\begin{eqnarray}\begin{array}{lclcl}
f(x) &= &\sum\limits_{\stackrel{n=-\infty}{n\mbox{ odd}}}^\infty \frac{2}{\pi i n}\{\rm\ e\}^{inx} &= &\frac{4}{\pi}\sum\limits_{\stackrel{n=1}{n\mbox{ odd}}}^\infty \frac{\{\rm\ e\}^{inx}-\{\rm\ e\}^{-inx}}{2in}\\
    &= &\frac{4}{\pi}\sum\limits_{\stackrel{n=1}{n\mbox{ odd}}}^\infty \frac{\sin (nx)}{n}. &&\end{array}
\end{eqnarray}
$$(ref-fourier-convergence-4)
```

## Explicit Summation of the First $N$ Terms

$$
\begin{eqnarray}
S_N(x) & = & \frac{4}{\pi}\sum\limits_{m=0}^N\frac{\sin
(2m+1)x}{(2m+1)}
   =\frac{4}{\pi} Im\left\{\sum\limits_{m=0}^N\frac{\{\rm\ e\}^{i(2m+1)x}}{(2m+1)}\right\}\\
S_N^\prime (x) & = & \frac{4}{\pi} Im\left\{\sum\limits_{m=0}^N i\{\rm\ e\}^{i(2m+1)x}\right\}\\
& = & \frac{4}{\pi} Im\left\{ i\{\rm\ e\}^{ix}\sum\limits_{m=0}^N {\big(\{\rm\ e\}^{i2x}\big)}^m\right\}\\
& = & \frac{4}{\pi} Im\left\{ i\{\rm\ e\}^{ix}\left(\frac{1+\{\rm\ e\}^{i2x}+\cdots +{\big(\{\rm\ e\}^{i2x}\big)}^N}{1-\{\rm\ e\}^{i2x}}\right) (1-\{\rm\ e\}^{i2x})\right\}\\
& = & \frac{4}{\pi} Im\left\{ i\{\rm\ e\}^{ix}\left(\frac{1-\{\rm\ e\}^{i2(N+1)x}}{1-\{\rm\ e\}^{i2x}}\right)\right\}\\
& = & \frac{4}{\pi}Im\left\{ i\left(\frac{1-\{\rm\ e\}^{i2(N+1)x}}{\{\rm\ e\}^{ix}-\{\rm\ e\}^{-ix}}\right)\right\}\\
& = & \frac{2}{\pi}Im\left\{\frac{\{\rm\ e\}^{i2(N+1)x}-1}{\sin x}\right\}\\
& = & \frac{2}{\pi}\frac{\sin 2(N+1)x}{\sin x}.
\end{eqnarray}
$$(ref-fourier-convergence-5)

Therefore

$$
\begin{eqnarray}\begin{array}{lclcl}
&&&&\hspace{.5in}\framebox{$t=2(N+1)u\qquad du=\frac{dt}{2(N+1)}$}\\
&&&&\hspace{.8in}\swarrow\\
S_N(x) &= &\frac{2}{\pi}\int\limits_0^x \frac{\sin 2(N+1)u}{\sin
u}\, du
   &\simeq &\frac{2}{\pi}\int\limits_0^{2(N+1)x} \frac{\sin t}{t}\, dt \end{array}
\end{eqnarray}
$$(ref-fourier-convergence-6)

````{prf:observation}
:label: observation-fourier-convergence-0
$\displaystyle S_N^\prime (x)=\frac{2}{\pi}\frac{\sin 2(N+1)x}{\sin
x}=0$ when $2(N+1)x_N=\pi$ thus the maximum value of $S_N(x)$ occurs
at

$$
\begin{equation}
\displaystyle x_N=\frac{\pi}{2(N+1)}
\end{equation}
$$(ref-fourier-convergence-7)

```{figure} ../img/fourier/half_range_derivative.png
:name: half_range_derivatives
:align: center

$(2/\pi) sin(2(N+1)x)/sin(x)$ for $N=5$
```

```{figure} ../img/fourier/half_range_gibbs.png
:name: half_range_gibbs
:align: center

Integral of $ (2/\pi) sin(2(N+1)x)/sin(x)$
```
````
