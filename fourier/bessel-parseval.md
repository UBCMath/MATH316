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

# Bessel's Inequality, Parseval's Theorem, Energy Convergence

In this lecture we consider the counterpart of Pythagoras' Theorem
for functions whose square is integrable. Square integrable
functions are associated with functions describing physical systems having finite energy. For a finite Fourier Series involving
$N$ terms we derive the so-called Bessel Inequality, in which
$N$ can be taken to infinity provided the function $f$ is square
integrable. The Bessel Inequality is shown to reduce to an equality
if and only if the Fourier Series $S_n(x)$ converges to $f$ in the
energy norm. The result is known as Parseval's Formula, which
has profound consequences for the completeness of the Fourier Basis
$\{1, \cos(\frac{n \pi x}{L}),\sin(\frac{n \pi x}{L})\}$. We see
that Parseval's Formula leads to a new class of sums for series of
reciprocal powers of $n$.

```{admonition} Key Concepts
Convergence of Fourier Series,
Bessel's Inequality, Paresval's Theorem, Plancherel theorem,
Pythagoras' Theorem, Energy of a function, Convergence in Energy,
completeness of the Fourier Basis.
```

## Bessel's Inequality

````{prf:definition}
Let $f(x)$ be a function that is
square-integrable on $[-L,L]$ i.e.,

$$
\begin{eqnarray*}
\int\limits_{-L}^L \big[ f(x)\big]^2\, dx<\infty ,
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-0)
in which case we write $f\in L_2[-L,L]$.
````

Consider the Fourier Series associated with $f(x)$, namely;

$$
\begin{eqnarray*}
f(x)\sim\frac{a_0}{2}+\sum\limits_{n=1}^\infty
a_n\cos\left(\frac{n\pi x}{L}\right) +b_n\sin\left(\frac{n\pi
x}{L}\right) =S_\infty
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-1)

Let

$$
\begin{eqnarray*}
S_N(x)=\frac{a_0}{2}+\sum\limits_{n=1}^N a_n\cos\left(\frac{n\pi
x}{L}\right) +b_n\sin\left(\frac{n\pi x}{L}\right) .
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-2)

Now

$$
\begin{eqnarray*}
\big[ f(x)-S_N(x)\big]^2=f^2(x)-2f(x)S_N(x)+S_N^2 (x)
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-3)

Consider the least-square error defined to be

$$
\begin{eqnarray*}
\mathcal{E}_2\big[ f,S_N\big] &=&\frac{1}{L}\int\limits_{-L}^L \big[
f(x)-S_N(x)\big]^2\, dx\\
&=&\frac{1}{L}\left\{\int\limits_{-L}^L f^2(x)\,
dx-2\int\limits_{-L}^L f(x)S_N(x)\, dx+\int\limits_{-L}^L S_N^2(x)\,
dx\right\}\\
&=&\frac{1}{L}\left\{ \langle f,f\rangle -2\langle f,S_N\rangle
+\langle S_N,S_N\rangle\right\}
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-4)

Now

$$
\begin{eqnarray*}
\langle S_N,S_N\rangle &=&\int\limits_{-L}^L
{\left[\frac{a_0}{2}+\sum\limits_{n=1}^N a_n\cos \left(\frac{n\pi
x}{L}\right) +b_n\sin\left(\frac{n\pi
x}{L}\right)\right]}^2\, dx\\
&=&\frac{a_0^2}{2}L+\sum\limits_{n=1}^N
a_n^2\int\limits_{-L}^L\cos^2 \left(\frac{n\pi x}{L}\right)\,
dx+b_n^2\int\limits_{-L}^L\sin^2 \left(\frac{n\pi x}{L}\right)\,
dx\\
&=& L\left[\frac{a_0^2}{2}+\sum\limits_{n=1}^N a_n^2+b_n^2\right]
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-5)

In addition,

$$
\begin{eqnarray*}
\langle f,S_N\rangle &=&\int\limits_{-L}^L f(x)S_N(x)\, dx\\
&=&\frac{a_0}{2}\int\limits_{-L}^L f(x)\, dx+\sum\limits_{n=1}^N
a_n\int\limits_{-L}^L f(x)\cos\left(\frac{n\pi x}{L}\right)\,
dx+b_n\int\limits_{-L}^L f(x)\sin\left(\frac{n\pi x}{L}\right)\,
dx\\
&=&\frac{a_0^2}{2}L+\sum\limits_{n=1}^N a_n^2L+b_n^2L.
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-6)

Therefore

$$
\begin{eqnarray*}
\mathcal{E}_2[f,S_N]=\frac{1}{L}\int\limits_{-L}^L {\big[
f(x)-S_N(x)\big]}^2\, dx=\frac{1}{L}\langle f,f\rangle
-\left\{\frac{a_0^2}{2}+\sum\limits_{n=1}^N a_n^2+b_n^2\right\}
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-7)

Now since $\displaystyle\mathcal{E}_2[f,S_N]=\int\limits_{-L}^L{\big[
f(x)-S_N(x)\big]}^2\, dx\geq 0$ it follows that

$$
\begin{eqnarray*}
\frac{a_0^2}{2}+\sum\limits_{n=1}^N
a_n^2+b_n^2\leq\frac{1}{L}\int\limits_{-L}^L f^2(x)\, dx=\frac{1}{L}
\langle f,f\rangle =E[f]
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-8)

where $E[f]$ is known as the energy of the $2L$-periodic function
$f$.

````{prf:theorem} Bessel's Inequality
Let $f\in L_2[-L,L]$ then

$$
\begin{eqnarray*}
\frac{a_0^2}{2}+\sum\limits_{n=1}^\infty
a_n^2+b_n^2\leq\frac{1}{L}\int\limits_{-L}^L f^2(x)\, dx
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-9)

in particular the series
$\displaystyle\frac{a_0}{2}+\sum\limits_{n=1}^\infty a_n^2+b_n^2$ is
convergent.
````

## Bessel's Inequality, Components of a Vector, and Pythagoras' Theorem

### 2D Analogue

Consider a 2D vector $f$, which is decomposed into components in
terms of two orthogonal unit vectors $\hat{e}_1$ and $\hat{e}_2$, i.e.

$$
\begin{eqnarray*}
\tilde{f}=a_1\hat{e}_1+a_2\hat{e}_2
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-10)

$$
\begin{eqnarray*}
|f|^2=\tilde{f}\cdot\tilde{f} &=&(a_1\hat{e}_1+a_2\hat{e}_2)\cdot
(a_1\hat{e}_1+a_2\hat{e}_2)\\
&=& a_1^2+a_2^2\mbox{ since $\hat{e}_k$ are orthogonal unit
vectors}\\
\mbox{Therefore }|f|^2 &=&a_1^2+a_2^2\mbox{ which is Pythagoras'
Theorem.}
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-11)

### 3D Analogue

Suppose we wish to expand a $3$-vector $\tilde{f}$
in terms of a set of $2$ basis vectors $\{\hat{e}_1,\hat{e}_2\}$.
Bessel's Inequality assumes the form

$$
\begin{eqnarray*}
a_1^2+a_2^2\leq |f|^2
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-12)

Since the subspace span $\{\hat{e}_1,\hat{e}_2\}$ (which represents
a plane in $\mathbb{R}^3$) does not include the whole of ${\mathbb
R}^3$ the vector $a_1\hat{e}_1+a_2\hat{e}_2\approx\tilde{f}$
represents the orthogonal projection of $\tilde{f}$ onto span
$\{\hat{e}_1,\hat{e}_2\}$. If we include the third basis vector
$\hat{e}_3$ in the basis, then the span
$\{\hat{e}_1,\hat{e}_2,\hat{e}_3\} ={\mathbb R}^3$. In this case the
set $\{\hat{e}_1,\hat{e}_2,\hat{e}_3\}$ are linearly independent and
of full rank and thus span the complete space ${\mathbb R}^3$.
$\{\hat{e}_1,\hat{e}_2,\hat{e}_3\}$ are in this case said to form a
complete set. In this case

$$
\begin{eqnarray*}
\tilde{f}=a_1\hat{e}_1+a_2\hat{e}_2+a_3\hat{e}_3
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-13)

and $\displaystyle |\tilde{f}|^2=a_1^2+a_2^2+a_3^2$ so that Bessel's
Inequality assumes the form of an equality, which in this trivial
case reduces to Pythagoras' Theorem. For a set of functions, that
are complete, the equivalent of Pythagoras' Theorem is Parseval's
Theorem.

## Parseval's Theorem

````{prf:theorem} Parseval's Identity
Let $f\in L_2[-L,L]$ then the Fourier coefficients $a_n$and $b_n$
satisfy Parseval's Formula

$$
\begin{eqnarray*}
\frac{a^2_0}{2}+\sum\limits_{n=1}^\infty
a_n^2+b_n^2=\frac{1}{L}\int\limits_{-L}^L f^2(x)\, dx=E[f]
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-14)

If and only if

$$
\begin{eqnarray*}
\lim\limits_{N\rightarrow\infty} \int\limits_{-L}^L\big[
f(x)-S_N(x)\big]^2\, dx=0.
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-15)
````

In this case the The Least Square Error assumes the form

$$
\begin{eqnarray}
\mathcal{E}_2[f,S_N] &=&\frac{1}{L}\int\limits_{-L}^L\big[
f(x)-S_N(x)\big]^2\, dx=\frac{1}{L}\int\limits_{-L}^L f^2(x)\,
dx-\left(\frac{a_0^2}{2}+\sum\limits_{n=1}^N a_n^2+b_n^2 \right)\nonumber\\
&=&\left(\frac{a_0}{2}+\sum\limits_{n=1}^\infty a_n^2+b_n^2\right)
-\left(\frac{a_0^2}{2}+\sum\limits_{n=1}^N a_n^2+b_n^2\right)\nonumber\\
&=&\sum\limits_{n=N+1}^\infty a_n^2+b_n^2
\end{eqnarray}
$$(ref-fourier-bessel-parseval-16)

### Parseval's Theorem for Odd Functions

````{prf:theorem} Parseval's Identity for Odd Functions

Let $\displaystyle f(x)=\sum\limits_{n=1}^\infty b_n \sin\left(\frac{n\pi
x}{L}\right)$ $0<x<L$. Then
$\displaystyle\frac{2}{L}\int\limits_0^L\big[ f(x)\big]^2\,
dx=\sum\limits_{n=1}^\infty b_n^2$.

**Proof:**

$$
\begin{eqnarray}
\int\limits_0^L\big[ f(x)\big]^2\, dx
   & = & \sum\limits_{m=1}^\infty \sum\limits_{n=1}^\infty b_m b_n\int\limits_0^L\sin
   \left(\frac{m\pi x}{L}\right)\sin\left(\frac{n\pi x}{L}\right)\, dx\\
& = & \sum\limits_{m=1}^\infty \sum\limits_{n=1}^\infty b_m b_n\cdot
\delta_{mn}\cdot\frac{L}{2}=\frac{L}{2}\sum\limits_{n=1}^\infty
b_n^2.
\end{eqnarray}
$$(ref-fourier-bessel-parseval-17)
````

````{prf:example}
:label: example-fourier-bessel-parseval-0
Recall for $x\in [0,2]$, $\displaystyle
f(x)=x=\frac{4}{\pi}\sum\limits_{n=1}^\infty\frac{(-1)^{n+1}}{n}\sin\left(\frac{n\pi
x}{2}\right)$. Therefore

$$
\begin{eqnarray}\begin{array}{lcrcl}
\frac{2}{L}\int\limits_0^L\big( f(x)\big)^2\, dx &= &\frac{2}{2}\int\limits_0^2 x^2\, dx &= &\left(\frac{4}{\pi}\right)^2 \sum\limits_{n=1}^\infty \frac{1}{n^2}\\
&\Rightarrow &\left.\frac{x^3}{3}\right|_0^2 &= &{\left(\frac{4}{\pi}\right)}^2\sum\limits_{n=1}^\infty \frac{1}{n^2}\\
&&\frac{\pi^2}{6} &= &\sum\limits_{n=1}^\infty
\frac{1}{n^2}\end{array}
\end{eqnarray}
$$(ref-fourier-bessel-parseval-18)
````

Note that $\displaystyle\sum\limits_{n=1}^\infty
\frac{1}{(2n)^2}=\frac{1}{2^2}\sum\limits_{n=1}^\infty\frac{1}{n^2}
=\frac{1}{4}\left(\frac{\pi^2}{6}\right)
=\frac{\pi^2}{24}$, and

$$
\begin{eqnarray}\begin{array}{rclcl}
&&\mbox{evens}&&\mbox{odds}\nonumber\\
\frac{\pi^2}{6}=\sum\limits_{n=1}^\infty \frac{1}{n^2}
   &= &\sum\limits_{m=1}^\infty \frac{1}{(2m)^2} &+ &\sum\limits_{m=0}^\infty\frac{1}{(2m+1)^2}\\
   &= &\frac{\pi^2}{24} &+ &\sum\limits_{m=0}^\infty \frac{1}{(2m+1)^2}\nonumber\end{array}
\end{eqnarray}
$$(ref-fourier-bessel-parseval-19)

therefore

$$
\begin{eqnarray}
\sum\limits_{m=0}^\infty \frac{1}{(2m+1)^2}=\frac{\pi^2}{6}-
\frac{\pi^2}{24}=\frac{\pi^2}{8}.
\end{eqnarray}
$$(ref-fourier-bessel-parseval-20)

```{note}
For Fourier Sine Components:

   $$
   \begin{eqnarray}
   \frac{2}{L}\int\limits_0^L \big( f(x)\big)^2\,
   dx=\sum\limits_{n=1}^\infty b_n^2.
   \end{eqnarray}
   $$(ref-fourier-bessel-parseval-21)
```

````{prf:example}
:label: example-fourier-bessel-parseval-1
Consider $f(x)=x^2$, $-\pi <x<\pi$.

The Fourier Series Expansion is:

$$
\begin{eqnarray}
x^2=\frac{\pi^2}{3}+4\sum\limits_{n=1}^\infty \frac{(-1)^n}{n^2}\cos
(nx).
\end{eqnarray}
$$(ref-fourier-bessel-parseval-22)

$$
\begin{eqnarray*}\begin{array}{rcccc}
n&1&\phantom{-}2&3&4\\ \cos\left(\frac{n\pi}{2}\right)
&0&-1&0&1\end{array}
\end{eqnarray*}
$$(ref-fourier-bessel-parseval-23)

Let

$$
\begin{eqnarray}\begin{array}{lclcl}
x=\frac{\pi}{2}&\Rightarrow &\frac{\pi^2}{4}&= &\frac{\pi^2}{3}+4\sum\limits_{n=1}^\infty \frac{(-1)^n}{n^2}\cos\left(\frac{n\pi}{2}\right)\\
&&-\frac{\pi^2}{12}&= &4\sum\limits_{k=1}^\infty
\frac{(-1)^k}{(2k)^2}\end{array}
\end{eqnarray}
$$(ref-fourier-bessel-parseval-24)

Therefore

$$
\begin{eqnarray}
\frac{\pi^2}{12} = \sum\limits_{k=1}^\infty \frac{(-1)^{k+1}}{k^2}.
\end{eqnarray}
$$(ref-fourier-bessel-parseval-25)

By Parseval's Formula:

$$
\begin{eqnarray}\begin{array}{rcl}
\frac{2}{\pi}\int\limits_0^\pi x^4\, dx &= &2{\left(\frac{\pi^2}{3}\right)}^2+16\sum\limits_{n=1}^\infty\frac{1}{n^4}\\
\frac{2}{\pi}\left.\frac{x^5}{5}\right|_0^\pi &=
&\frac{2\pi^4}{9}+16\sum\limits_{n=1}^\infty
\frac{1}{n^4}\end{array}\qquad\begin{array}{c}\frac{9-5}{45}=\frac{4}{45}=\frac{8}{90}\\
\frac{1}{90}\end{array}\end{eqnarray}
$$(ref-fourier-bessel-parseval-26) 

Therefore

$$
\begin{eqnarray}
\frac{\pi^4}{90} &= &\sum\limits_{n=1}^\infty\frac{1}{n^4}=
\zeta(4),
\end{eqnarray}
$$(ref-fourier-bessel-parseval-27)

where $\zeta$ is the Riemann Zeta Function defined by:

$$
\begin{eqnarray}  \zeta(s)&=&\sum\limits_{n=1}^\infty\frac{1}{n^s}, \
s=\sigma+\rm(i) \tau, \ \sigma =\rm{Re\{ s\}}>1
\end{eqnarray}
$$(ref-fourier-bessel-parseval-28)
````
