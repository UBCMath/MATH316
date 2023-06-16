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

# Bessel's Equation and Bessel Functions

In this lecture we will consider the Frobenius series solution of
the Bessel equation, which arises during the process of separation
of variables for problems with radial or cylindrical symmetry.
Depending on the parameter $\nu$ in  Bessel's equation, we obtain
roots of the indicial equation that are: distinct and real,
repeated, and which differ by an integer.

````{admonition} Key Concepts
Frobenius Series Solutions, Bessel's equation; Bessel Functions.
````

## Bessel Functions

### Bessel's Function of Order $\nu\notin\{ \ldots,-2,-1,0,1,2 \ldots\}$

$$
\begin{equation}
Ly=x^2y^{\prime\prime}+xy^\prime +(x^2-\nu^2)y=0
\end{equation}
$$(ref-odes-bessel-0)

$x=0$ is a regular Singular Point: therefore let $\displaystyle
y=\sum\limits_{n=0}^\infty a_nx^{n+r}$.

$$
\begin{eqnarray}
0=\sum\limits_{n=0}^\infty a_n\big[ (n+r)(n+r-1)+(n+r)-\nu^2\big] x^{n+r} &
+ &\sum\limits_{n=0}^\infty a_nx^{n+r+2} \\
& &\, m=n+2\quad n=m-2  \notag \\
& &\, n=0\Rightarrow m=2  \notag
\end{eqnarray}
$$(ref-odes-bessel-1)

$$
\begin{eqnarray}
0=\sum\limits_{m=2}^\infty\left\{ a_m\big[ (m+r)^2-\nu^2\big]
+a_{m-2}\right\} x^{m+r} & + & a_0\left\{ r^2-\nu^2\right\} x^r \\
& &\quad +\, a_1\left\{ (1+r)^2-\nu^2\right\} x^{r+1}  \notag
\end{eqnarray}
$$(ref-odes-bessel-2)

$$
\begin{eqnarray}
\begin{array}{ll}
x^r> & a_0\not= 0\Rightarrow r=\pm\nu\quad\mbox{Indicial Eq. Roots} \\
x^{r+1}> & a_1\big\{ r^2+2r+1-\nu^2\big\} =0,\quad
\mbox{$a_1=0$ provided $\nu\not=-\frac{1}{2}$ and if $\nu=-\frac{1}{2}$ then $a_1$ is arbitrary.} \\
x^{m+r}> & a_m=\displaystyle -\frac{a_{m-2}}{(m+r)^2-\nu^2}\quad m\geq 2
\end{array}
\end{eqnarray}
$$(ref-odes-bessel-3)

$r=\nu$:

$$
\begin{eqnarray}
a_m & = &-\frac{a_{m-2}}{(m+\nu )^2-\nu^2}=-\frac{a_{m-2}}{m^2+2m\nu}=-\frac{
a_{m-2}}{m(m+2\nu )}  \notag \\
a_2=-\frac{a_0}{2(2+2\nu )} & = & -\frac{a_0}{2^2(1+\nu )}\quad a_4=-\frac{
a_2}{4(4+2\nu )}=\quad \frac{(-1)^2a_0}{2.2^4(2+\nu )(1+\nu )}  \notag \\
& & \quad \ldots a_{2m}=\frac{(-1)^ma_0}{m!2^{2m}(1+\nu )\ldots (m+\nu )} \\
y_1(x) & = & x^\nu \sum\limits_{m=0}^\infty \frac{(-1)^m(x/2)^{2m}}{m!(1+\nu
)(2+\nu )\ldots (m+\nu )}\quad \overset{x\rightarrow 0}{\rightarrow 0}
\notag
\end{eqnarray}
$$(ref-odes-bessel-4)

$r=-\nu$:

$$
\begin{eqnarray}
a_m & = & -\frac{a_{m-2}}{m(m-2\nu )}  \notag \\
a_2 = -\frac{a_0}{2(2-2\nu )}& = & -\frac{a_0}{2^2(1-\nu )},\quad a_4=-\frac{
a_2}{4(4-2\nu )}=\frac{(-1)^2a_0}{22^4(1-\nu )(2-\nu )}  \notag \\
& &\quad \ldots a_{2m}=\frac{(-1)^ma_0}{m!2^{2m}(1-\nu )\ldots (m-\nu )} \\
y_2(x) & = & x^{-\nu}\sum\limits_{m=0}^\infty\frac{(-1)^m(x/2)^{2m}}{
m!(1-\nu )\ldots (m-\nu )}\quad\overset{x\rightarrow 0}{\rightarrow\infty}
\notag
\end{eqnarray}
$$(ref-odes-bessel-5)

### Bessel's Function of Order $\nu =0$ - Repeated Roots

In this case

$$
Ly=x^2y+xy^\prime +x^2y=0
$$

$$
\begin{eqnarray}
y & = & \sum\limits_{n=0}^\infty a_nx^{n+r}  \notag \\
Ly & = & \sum\limits_{n=0}^\infty a_n\big\{ (n+r)(n+r-1)+(n+r)\big\}
x^{n+r}+a_nx^{n+r+2}=0  \notag \\
& &\hspace{3in}m=n+2\quad n=m-2 \\
0 & = & \sum\limits_{n=2}^\infty \big[ a_n(n+r)^2+a_{n-2}\big] x^{n+r}+a_0
\big[ r(r-1)+r\big] x^r+a_1\big[ (r+1)r+r+1\big] x^{r+1}=0  \notag
\end{eqnarray}
$$(ref-odes-bessel-6)

The indicial equation is: $a_0r^2=0\quad r_{1,2}=0,0\quad$ a double root.

$r_1=0\Rightarrow a_1.1=0\Rightarrow a_1=0$.

__Recursion:__ $a_n=\displaystyle -\frac{a_{n-2}}{(n+r)^2}\quad n\geq 2$.

$$
\begin{eqnarray}
\begin{array}{c}
a_2=-\displaystyle\frac{a_0}{2^2}; \quad a_4=-\frac{a_2}{4^2}=\frac{a_0}{
2^24^2};\quad a_6=-\frac{a_4}{6^2}=-\frac{a_0}{2^24^26^2};\quad a_8=\frac{a_0
}{2^24^26^28^2}
\end{array}
\end{eqnarray}
$$(ref-odes-bessel-7)

$$
\begin{eqnarray}
a_{2m} & = & \frac{(-1)^m}{2^{2m}(m!)^2}a_0 \\
y_1(x) & = & \left\{ 1+\sum\limits_{m=1}^\infty\frac{(-1)^mx^{2m}}{
2^{2m}(m!)^2}\right\} = J_0(x)  \notag
\end{eqnarray}
$$(ref-odes-bessel-8)

```{figure} ../img/odes/zeroth_order_bessel.png
:name: bessel-functions
:alt: bessel-functions
:align: center

Zeroth order bessel functions $j_0(x)$ and $Y_0(x)$
```

To get a second solution

$$
\begin{eqnarray}
y(x,r) & = & a_0x^r\left\{ 1-\frac{x^2}{(2+r)^2}+\frac{x^4}{(2+r)^2(4+r)^2}
+\cdots +\frac{(-1)^mx^{2m}}{(2+r)^2(4+r)^2\ldots (2m+r)^2}\right.  \notag \\
& &\quad +\, \phantom{\bigg\{} \cdots \bigg\} \\
\left. \frac{\partial y}{\partial r}(x,r)\right| _{r=r_1} & = & a_0\log
xy_1(x)+a_0x^r\sum\limits_{m=1}^\infty (-1)^mx^{2m}\frac{\partial}{\partial r
}\left\{\frac{1}{(2+r)^2\ldots (2m+r)^2}\right\} .  \notag
\end{eqnarray}
$$(ref-odes-bessel-9)

Let

$$
\begin{eqnarray}
a_{2m}(r)=\left\{ \,\,\,\right\}\Rightarrow \ln a_{2m}(r) & = & -2\ln (2+r)
- \ldots -2\ln (2m+r) \\
a_{2m}^\prime (0) & = & \left. \left( -\frac{2}{2+r}-\frac{2}{4+r}\cdots -
\frac{2}{(2m+r)}\right)\right|_{r=0} a_{2m}(0)  \notag \\
& = & \left( -1-\frac{1}{2}- \ldots -\frac{1}{m}\right)
a_{2m}(0)=-H_ma_{2m}(0).  \notag
\end{eqnarray}
$$(ref-odes-bessel-10)

Let $H_m=1+\displaystyle\frac{1}{2}+\cdots +\frac{1}{m}$. Therefore

$$
\begin{equation}
y_2(x)=J_0(x)\ln x+\sum\limits_{m=1}^\infty \frac{(-1)^{m+1}H_m}{2^{2m}(m!)^2
}x^{2m}\quad x>0.
\end{equation}
$$(ref-odes-bessel-11)

It is conventional to define

$$
\begin{equation}
Y_0(x)=\frac{2}{\pi}\big[ y_2(x)+(\gamma -\log 2)J_0(x)\big]
\end{equation}
$$(ref-odes-bessel-12)

where

$$
\begin{eqnarray}
\gamma & = & \lim\limits_{n\rightarrow\infty} (H_n-\log n)=0{.}5772\quad
\mbox{Euler's Constant}  \notag \\
y(x) & = & c_1J_0(x)+c_2Y_0(x).
\end{eqnarray}
$$(ref-odes-bessel-13)

### Bessel's Function of Order $\nu =\frac{1}{2}$

Consider the case $\nu=1/2$: $Ly=x^2y^{\prime\prime} +xy^\prime +\displaystyle
\left( x^2-\frac{1}{4}\right) \quad y=0$. Let

$$
\begin{equation}
y=\sum\limits_{n=0}^\infty a_nx^{n+r}
\end{equation}
$$(ref-odes-bessel-14)

$$
\begin{eqnarray}
Ly & = & \sum\limits_{n=0}^\infty a_n\left\{ (n+r)^2-\frac{1}{4}\right\}
x^{n+r}+\sum\limits_{n=0}^\infty a_nx^{n+r+2}=0
\begin{array}{lcl}
m & = & n+2 \\
n & = & m-2 \\
n & = & 0\Rightarrow m=2
\end{array}
\\
Ly & = & a_0\left\{ r^2-\frac{1}{4}\right\} x^r +a_1\left\{ (r+1)^2-\frac{1}{4}
\right\} x^{r+1} +\sum\limits_{n=2}^\infty \left[ a_n\left\{ (n+r)^2-\frac{1}{4}
\right\} +a_{n-2}\right] x^{n+r}=0.  \notag
\end{eqnarray}
$$(ref-odes-bessel-15)

__Indicial Equation:__ $r^2-\displaystyle\frac{1}{4}=0,\quad r=\pm
\displaystyle\frac{1}{2}\quad$ Roots differ by an integer.

__Recurrence:__ $a_n=-\displaystyle\frac{a_{n-2}}{(n+r)^2-\frac{1}{4}}
\quad n\geq 2$.

$r_1=+ 1/2$:

$$
\begin{eqnarray}
\begin{array}{c}
a_n=-\displaystyle\frac{a_{n-2}}{(n+\frac{1}{2})^2-\frac{1}{4}}=-\frac{
a_{n-2}}{(n+1)n}\quad n\geq 2;\left(\frac{9}{4}-\frac{1}{4}\right)
a_1=0\Rightarrow a_1=0 \\
a_2=-\displaystyle\frac{a_0}{3.2}\quad a_4=\frac{(-1)^2a_0}{5.4.3.2}\ldots
a_{2n}=\frac{(-1)^na_0}{(2n+1)!} \\
y_1(x)=\displaystyle x^{\frac{1}{2}}\sum\limits_{n=0}^\infty \frac{
(-1)^nx^{2n}}{(2n+1)!}=x^{-\frac{1}{2}}\sum\limits_{n=0}^\infty \frac{
(-1)^nx^{2n+1}}{(2n+1)!}=x^{-\frac{1}{2}}\sin x
\end{array}
\end{eqnarray}
$$(ref-odes-bessel-16)

$r_2=-\displaystyle\frac{1}{2}$:

$$
\begin{eqnarray}
\begin{array}{c}
a_n=-\displaystyle\frac{a_{n-2}}{(n-\frac{1}{2})^2-\frac{1}{4}}= - \frac{
a_{n-2}}{n(n-1)},\quad n\geq 2, \\
n=1\Rightarrow a_1\displaystyle\left\{\left( -\frac{1}{2}+1\right)^2 -\frac{1
}{4}\right\} = a_1.0=0\quad a_1\mbox{ and }a_0\mbox{
arbitrary}.
\end{array}
\end{eqnarray}
$$(ref-odes-bessel-17)

$a_0$:

$$
\begin{equation}
a_2= - \frac{a_0}{2.1}\quad a_4=\frac{(-1)^2a_0}{4.3.2.1}\ldots\quad a_{2n}=
\frac{(-1)^na_0}{(2n)!}
\end{equation}
$$(ref-odes-bessel-18)

$a_1$:

$$
\begin{equation}
a_3= - \frac{a_1}{3.2}\quad a_5=\frac{(-1)^2a_1}{5.4.3.2}\quad a_{2n+1}=
\frac{(-1)^na_1}{(2n+1)!}
\end{equation}
$$(ref-odes-bessel-19)

$$
\begin{eqnarray}
y_2(x) & = & a_0x^{-\frac{1}{2}}\sum\limits_{n=0}^\infty \frac{(-1)^nx^{2n}}{
(2n)!} +a_1x^{-\frac{1}{2}}\sum\limits_{n=0}^\infty \frac{(-1)^nx^{2n+1}}{
(2n+1)!}  \notag \\
& = & a_0x^{-\frac{1}{2}}\cos x+a_1x^{-\frac{1}{2}}\sin x \\
& &\hspace{1.2in}\nwarrow\mbox{ included in $y_1(x)$}.  \notag
\end{eqnarray}
$$(ref-odes-bessel-20)

````{note}
In this case the recursion spawns another solution
for the smaller root $r=-\frac{1}{2}$ so we get away without having
to do anything special to get another solution. In the next
subsection we give an example where this is not the case and we have
to use our differentiation with respect to $r$ trick. We could
always use the method of reduction of order along with the first
solution.
````

### The Roots Differ by an Integer - An Example for Enrichment

````{prf:example}
:label: example-odes-bessel-0
Let $Ly=xy^{\prime\prime}-y=0,\quad x=0$ is a regular singular point.

$$
\begin{eqnarray}
& &y=\sum\limits_{n=0}^\infty c_nx^{n+\alpha}  \notag \\
& &\sum\limits_{n=0}^\infty c_n(n+\alpha )(n+\alpha -1)x^{n+\alpha
-1}-c_nx^{n+\alpha}=0  \notag \\
& &\hspace{2.7in}\uparrow \\
& &\hspace{2.6in}p-1=n  \notag \\
& &\sum\limits_{n=1}^\infty \left\{ c_n(x+\alpha )(n+\alpha
-1)-c_{n-1}\right\} x^{n+\alpha -1}+c_0(\alpha -1)\alpha x^{\alpha -1}=0
\notag
\end{eqnarray}
$$(ref-odes-bessel-21)

__Indicial Equation:__ $(\alpha -1) \alpha =0 \ \Rightarrow \
\alpha =0,1\quad$ differ by integer.

__Recurrence Relation:__ $c_n=\displaystyle\frac{c_{n-1}}{(n+\alpha
)(n+\alpha -1)}\quad n\geq 1$.

```{note}
When $\alpha = 0$, $c_1$ blows up!
```

Let $\alpha =1\Rightarrow\quad c_1=\displaystyle\frac{c_0}{2},\quad c_2=
\frac{c_0}{12},\ldots $.

$$
\begin{equation}
y_1(x)=c_0x \left( 1+\frac{x}{2}+\frac{x^2}{12}+\cdots \right) =c_0u_1(x).
\end{equation}
$$(ref-odes-bessel-22)

__Second Solution:__

$$
\begin{eqnarray}
\bar{y}(x,\alpha ) & = & \alpha y(x,\alpha )=c_0x^\alpha\left\{\alpha +\frac{x}{
1+\alpha}+\frac{x^2}{(1+\alpha )(2+\alpha )(1+\alpha )}+\cdots\right\}
\notag \\
\frac{\partial\bar{y}}{\partial\alpha} & = & c_0x^\alpha\ln x\left\{\alpha +
\frac{x}{1+\alpha}+\cdots\right\} \\
& &\quad +\, c_0x^\alpha \left\{ 1-\frac{x}{(1+\alpha )^2}-\frac{x^2}{(1+\alpha
)^2(2+\alpha )} \left[\frac{2}{(1+\alpha )}+\frac{1}{(2+\alpha )}\right]
+\cdots\right\}  \notag \\
\left. \frac{\partial\bar{y}}{\partial\alpha}\right|_{\alpha =0}& = &
c_0\left\{ x+\frac{x^2}{2}+\frac{x^3}{12}+\cdots\right\}\ln x+c_0 \left\{
1-x-\frac{5}{4}x^2-\cdots\right\} = c_0u_2.  \notag \\
\mbox{Therefore }y(x) & = & (A+B\ln x)\left( x+\frac{x^2}{2}+\frac{x^3}{12}
+\cdots \right) + B\left( 1-x-\frac{5}{4}x^2-\cdots\right) .  \notag
\end{eqnarray}
$$(ref-odes-bessel-23)
````
