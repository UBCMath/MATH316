
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

# Sturm-Liouville Problems Involving The Cauchy-Euler Equation - Applications

In this lecture we look at eigenvalue problems involving
equidimensional or Cauchy-Euler differential operators.

```{admonition} Key Concepts
Eigenvalue Problems, Sturm-Liouville
Boundary Value Problems; Cauchy-Euler Equations; Equidimensional
equations.
```

Reference Section: Boyce and Di Prima Section 11.1 and 11.2

## Variable Coefficient Bvp - Eigenfunctions Involving Solutions To The Euler Equation:

````{prf:example}
:label: example-sturm-liouville-cauchy-euler-0 Eigenfunctions Involving Solutions to an Euler Equation

$$
\begin{eqnarray}
\begin{array}{ll}
{(x^2\phi^\prime )}^\prime +\lambda\phi =0\quad &1<x<2\\
\phi (1)=0, \, \, \phi (2)=0&\\
x^2\phi^{\prime\prime}+2x\phi^\prime +\lambda\phi =0 &\mbox{A
Cauchy-Euler Eq.}\end{array}
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-0)

Let

$$
\begin{eqnarray}
\phi (x)=x^r\quad r(r-1)+2r+\lambda =r^2+r+\lambda =0.
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-1)

Therefore

$$
\begin{eqnarray}
r=\frac{-1\pm\sqrt{1-4\lambda}}{2}=r_1,r_2.
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-2)

**Case I - $\mathbf{\lambda =\frac{1}{4}}$:**

$$
\begin{eqnarray}
\phi (x) &=& c_1x^{-\frac{1}{2}}+c_2x^{-\frac{1}{2}}\log x\\
\phi (1) &=& c_1=0\,\, \phi (2)=c_22^{-\frac{1}{2}}\log 2
=0\Rightarrow c_2=0
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-3)
so there is no nontrivial eigenfunction for $\lambda =1/4$.

**Case II - $\mathbf{\lambda \ne\frac{1}{4}}$:**

$$
\begin{eqnarray}
\phi (x) &=& c_1x^{r_1}+c_2x^{r_2}\\
\phi (1)&=&c_1+c_2=0\,\, c_2=-c_1\\
\phi (2)&=& c_1\big( 2^{r_1}-2^{r_2}\big) =0
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-4)

$$
\begin{eqnarray}\begin{array}{l}
2^{r_1-r_2}=1\\
 e^{(r_1-r_2)\ln 2}=1= e^{2\pi in}\\
(r_1-r_2)\ln 2=2\pi in\\
r_1-r_2=\sqrt{1-4\lambda}=2\pi ni/\ln (2)\end{array}
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-5)

Thus to obtain nontrivial solutions we require $1-4\lambda <0$ which
implies \framebox{$\lambda>\frac{1}{4}$}. Thus for $\lambda >\frac{1}{4}$

$$
\begin{eqnarray}
\sqrt{1-4\lambda}=i\sqrt{4\lambda -1}=2\pi ni/\ln (2).
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-6)

The Eigenvalues are:

$$
\begin{eqnarray}
\lambda_n =\frac{1}{4}+\frac{\pi^2 n^2}{[\ln (2)]^2},\, 4\lambda_n
-1=\frac{4\pi^2n^2}{[\ln (2)]^2}=(2\beta_n)^2\quad \beta_n=(n\pi
/\ln 2).
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-7)

The corresponding roots $r_1$ and $r_2$ are as follows

$$
\begin{eqnarray}
{(r_1)}_n &=& - \frac{1}{2}+i\beta_n \mbox{ and }{(r_2)}_n=- \frac{1}{2} -i\beta_n\\
\phi_n(x) &=& c_nx^{-\frac{1}{2}}\left( x^{i\beta_n} - x^{-i\beta_n}\right)\\
&=& c_nx^{ -\frac{1}{2}}\left[ e^{i\beta_n \ln x}
- e^{-i\beta_n\ln x}\right]\\
&=& d_nx^{ -\frac{1}{2}}\sin \left(\beta_n\ln x\right)\\
&=& d_nx^{-\frac{1}{2}}\sin \left[ n\pi\frac{\ln x}{\ln (2)}\right]
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-8)

Now let us consider expanding a function $f(x)$ in terms of a
'Fourier Series' of these new eigenfunctions in the following form

$$
\begin{eqnarray}
f(x)=\sum\limits_{n=1}^{\infty }c_{n}\phi
_{n}(x)=\sum\limits_{n=1}^{\infty }c_{n}x^{-1/2}\sin \left(
\frac{n\pi \ln x}{\ln 2}\right) \label{eq_genFS}
\end{eqnarray}
$$(eq_genFS)

In order to determine the coefficients $c_n$ we project the
function $f(x)$ onto the basis functions $\phi_n(x)$ as follows:

$$
\int\limits_{1}^{2}f(x)x^{-1/2}\sin \left( \frac{m\pi \ln x}{\ln
2}\right) dx=\sum\limits_{n=1}^{\infty
}c_{n}\int\limits_{1}^{2}x^{-1/2}\sin \left(
\frac{m\pi \ln x}{\ln 2}\right) x^{-1/2}\sin \left( \frac{n\pi \ln x}{\ln 2}
\right) dx
$$

The integrals under the summation can be evaluated by making
the simple substitution $z=\ln x$ \ so that $dz=\frac{dx}{x}$, in
which case:

$$
\int\limits_{1}^{2}\sin \left( \frac{m\pi \ln x}{\ln 2}\right) \sin
\left( \frac{n\pi \ln x}{\ln 2}\right)
\frac{dx}{x}=\int\limits_{0}^{\ln 2}\sin \left( \frac{m\pi z}{\ln
2}\right) \sin \left( \frac{n\pi z}{\ln 2}\right) dz=\delta
_{mn}\frac{\ln 2}{2}
$$

Substituting this result into {eq}`eq_genFS` we obtain the
following expression for the Fourier coefficients $c_n$:

$$
c_{n}=\frac{2}{\ln 2}\int\limits_{1}^{2}f(x)x^{-1/2}\sin \left(
\frac{n\pi \ln x}{\ln 2}\right) dx
$$
````

````{prf:example}
:label: example-sturm-liouville-cauchy-euler-1 A Variable Coefficient Heat Conduction Problem with Cauchy-Euler Eigenfunctions

$$
\begin{eqnarray}\begin{array}{ll}
u_t=D{(x^2u_x)}_x-u\quad &1<x<2\quad t>0\\
u(1,t)=0=u(2;t)&u(x,0)=f(x).\end{array}
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-10)

Let

$$
\begin{eqnarray}
u(x,t)&=&X(x)T(t)\\
\frac{\dot{T}(t)}{DT(t)}&=&\frac{{(x^2X^\prime )}^\prime}{X}-\frac{1}{D}\\
\frac{\dot{T}(t)}{DT(t)}+\frac{1}{D}&=&\frac{{(x^2X^\prime )}^\prime}{X}=-\lambda
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-11)

$$
\begin{eqnarray}\begin{array}{ll}
\dot{T}+(1+D\lambda )T=0 \quad &T(t)=c e^{-(1+D\lambda )t}\\
\left.\begin{array}{c}{(x^2X^\prime )}^\prime +\lambda X=0\\
X(1)=0=X(2)\end{array}\right\} \quad &\begin{array}{l}\lambda_n =\frac{1}{4}+\frac{(\pi n)^2}{[\ln (2)]^2}; X_n(x)=x^{-\frac{1}{2}}\sin\left( n\pi\frac{\ln x}{\ln 2}\right)\\
n=1,2,\ldots \end{array}\end{array}\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-12)

$$
\begin{eqnarray}
u(x,t)&=& x^{-\frac{1}{2}}\sum\limits_{n=1}^\infty c_n e^{-(1+D\lambda_n)t}\sin\left( n\pi\frac{\ln x}{\ln 2}\right)\\
f(x)&=& u(x,0)=x^{-\frac{1}{2}}\sum\limits_{n=1}^\infty
c_n\sin\left( n\pi\frac{\ln x}{\ln 2}\right)
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-13)

where

$$
c_{n}=\frac{2}{\ln 2}\int\limits_{1}^{2}f(x)x^{-1/2}\sin \left(
\frac{n\pi \ln x}{\ln 2}\right) dx
$$
````

### Solving Laplace'S Equation Using Cauchy-Euler Eigenfunctions

```{figure} ../img/sturm-liouville/annular.png
:name: annular
:align: center

Annular sector subtending an arc of $\alpha$ radians
between the radii $a$ and $b$: $1=b<a=2$.
```

$$
\begin{eqnarray}
\Delta u&=& u_{rr}+\frac{1}{r}u_r+\frac{1}{r^2}u_{\theta\theta}=0\quad 1<r<2,\  0<\theta<\alpha\phantom{\alpha\alpha}\\
u(r,0)&=&0,\ u(r,\alpha )=f( r )\\
u(1,\theta )&=&0, u(2,\theta )=0\\
u(r,\theta )&=&R( r )\Theta (\theta )\\
\frac{r^2R^{\prime\prime}+rR^\prime}{R( r )}&=&
-\frac{\Theta^{\prime\prime}}{\Theta }=-\lambda^2\quad \mbox{(because of Homog. BC)}
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-14)

$\mathbf{\Theta }$: $\displaystyle\begin{array}{ll}
\Theta^{\prime\prime}-\lambda^2\Theta =0\quad &\Theta =C\cosh\lambda\theta +D\sinh\lambda\theta\\
\Theta (0)=0 &\Theta (0)=C=0\Rightarrow\Theta (\theta
)=D\sinh\lambda\theta\end{array}$

$\mathbf{R}$: $r^2R^{\prime\prime}+rR^\prime +\lambda^2 R=0
(\star), \quad R(1)=0=R(2)$ Although we can easily see that dividing
through by $r$ we can reduce $(\star )$ to S-L form, let us use the
integrating factor

$$
\begin{eqnarray}
\mu(r) =\frac{1}{P} e^{\int\frac{Q}{P}dr}=\frac{1}{r^2}
e^{\int\frac{r}{r^2}dr}=\frac{1}{r^2} e^{\ln r}=\frac{1}{r}.
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-15)

Therefore $-\frac{1}{r}\cdot (\star )$ and some rearrangement
implies:

$$
-rR^{\prime \prime }-R^{\prime }=-\left( rR^{\prime }\right) ^{\prime }=
\frac{\lambda ^{2}}{r}R
$$

Now let us look for Eigenvalues and Eigenfunctions to $(\star )$.
Let

$$
\begin{eqnarray}
R( r ) &=& r^\gamma\Rightarrow\gamma (\gamma -1)+\gamma +\lambda^2 =\gamma^2 +\lambda ^2=0\quad\gamma =\pm i\lambda .\\
R( r ) &=& c_1r^{i\lambda }+c_2r^{-i\lambda}\quad r^{i\lambda }= e^{i\lambda \ln r}\\
  &=& A\cos (\lambda\ln r)+B\sin (\lambda \ln r)\\
R(1) &=& A\cos\left[\lambda (\ln 1)\right] +B\sin (\lambda\ln 1)=A=0\\
R(2) &=& B\sin \left[\lambda\ln 2\right] =0\Rightarrow\lambda_n \ln
2=n\pi\quad n=1,2,\ldots
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-16)

and the corresponding Eigenfunctions are $R_n=\sin\left(
n\pi\frac{\ln r}{\ln 2}\right) $. Therefore

$$
\begin{eqnarray}
u(r,\theta ) =\sum\limits_{n=1}^\infty B_n\sinh(\lambda _n\theta )\sin
(\lambda_n \ln r).
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-17)
Now match BC:

$$
\begin{eqnarray}
f( r )=u(r,\alpha )=\sum\limits_{n=1}^\infty B_n\sinh (\lambda_n \alpha
)\sin\left( n\pi\frac{\ln r}{\ln 2}\right)\,
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-18)

Now making the substitution $x=\ln r$ so that $dx=\frac{dr}{r}$ we
can reduce the orthogonality integrals for the Fourier coefficients
to the form:

$$
\begin{eqnarray}
\int\limits_1^2\frac{1}{r}\sin\left(\frac{m\pi\ln r}{\ln
2}\right)\sin\left(\frac{n\pi\ln r}{\ln 2}\right)\, dr
=\int\limits_0^{\ln 2}\sin\left(\frac{m\pi x}{\ln
2}\right)\sin\left(\frac{n\pi x}{\ln 2}\right)\, dx
=\left\{\begin{array}{ll}0\quad &m\ne n\\ \frac{\ln
2}{2}&m=n\end{array}\right. .
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-19)
Therefore

$$
\begin{eqnarray}
B_n=\frac{2}{\ln 2\sinh\left(\frac{n\pi\alpha}{\ln
2}\right)}\int\limits_1^2 \frac{f( r )}{r}\sin\left(\frac{n\pi \ln
r}{r}\right)\, dr.
\end{eqnarray}
$$(ref-sturm-liouville-cauchy-euler-20)
