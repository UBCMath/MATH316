
$$
\newenvironment{remark}{\begin{Remark}\rm}{\end{Remark}}
\newenvironment{exercise}{\begin{Exercise}\rm}{\end{Exercise}}
\newcommand{\N}[1]{\left\|#1\right\|}
\newcommand{\abs}[1]{|#1|}
\newcommand{\mat}[1]{{\mathbf #1}}
\newcommand{\vect}[1]{\underline{#1}}
\newcommand{\njump}[1]{[|#1|]}
\newcommand{\pa}{\partial}
\newcommand{\De}{\Delta}
\newcommand{\ra}{\rightarrow}
\newcommand{\dst}{\displaystyle}
\newcommand{\la}{\lambda}
\newcommand{\al}{\alpha}
\newcommand{\rme}{{\rm e}}
\newcommand{\rmi}{{\rm i}}

\newcommand{\R}{\mathbb{R}}
\newcommand{\bke}[1]{\left ( #1 \right )}
\newcommand{\bkt}[1]{\left [ #1 \right ]}
\newcommand{\bket}[1]{\left \{ #1 \right \}}
\newcommand{\norm}[1]{\left \| #1 \right \|}
\newcommand{\bka}[1]{\left \langle #1 \right \rangle}
\renewcommand{\th}{\theta}
\newcommand{\om}{\omega}
\newcommand{\pd}{\partial}
\newcommand{\dis}{\displaystyle}
\newcommand{\ve}[1]{\mathbf{#1}}
\newcommand{\what}[1]{\widehat{#1}}
\newcommand{\tint}{\int\kern-.6em\int\kern-.6em\int}
$$

# Sturm-Liouville Problems Involving The Cauchy-Euler Equation - Applications

In this lecture we look at eigenvalue problems involving
equidimensional or Cauchy-Euler differential operators.

```{admonition} Key Concepts
Eigenvalue Problems, Sturm-Liouville
Boundary Value Problems; Cauchy-Euler Equations; Equidimensional
equations.
```

\vskip 0.1 in
 Reference Section: Boyce and Di Prima Section 11.1 and 11.2

## Variable Coefficient Bvp - Eigenfunctions Involving Solutions To The Euler Equation:
\begin{example}  Eigenfunctions involving solutions to an Euler Equation:

$$
\begin{eqnarray}
\begin{array}{ll}
{(x^2\phi^\prime )}^\prime +\la\phi =0\quad &1<x<2\\
\phi (1)=0, \, \, \phi (2)=0&\\
x^2\phi^{\prime\prime}+2x\phi^\prime +\la\phi =0 &\mbox{A
Cauchy-Euler Eq.}\end{array}
\end{eqnarray}
$$(ref0)
\end{example}
Let

$$
\begin{eqnarray}
\phi (x)=x^r\quad r(r-1)+2r+\la =r^2+r+\la =0.
\end{eqnarray}
$$(ref1)
Therefore

$$
\begin{eqnarray}
r=\frac{-1\pm\sqrt{1-4\la}}{2}=r_1,r_2.
\end{eqnarray}
$$(ref2)


{Case I:} \framebox{$\mathbf{\la =\frac{1}{4}}$:}

$$
\begin{eqnarray}
\phi (x) &=& c_1x^{-\frac{1}{2}}+c_2x^{-\frac{1}{2}}\log x\\
\phi (1) &=& c_1=0\,\, \phi (2)=c_22^{-\frac{1}{2}}\log 2
=0\Rightarrow c_2=0
\end{eqnarray}
$$(ref3)
so there is no nontrivial eigenfunction for $\la =1/4$.

{Case II:} \framebox{$\mathbf{\la \ne\frac{1}{4}}$:}

$$
\begin{eqnarray}
\phi (x) &=& c_1x^{r_1}+c_2x^{r_2}\\
\phi (1)&=&c_1+c_2=0\,\, c_2=-c_1\\
\phi (2)&=& c_1\big( 2^{r_1}-2^{r_2}\big) =0
\end{eqnarray}
$$(ref4)

$$
\begin{eqnarray}\begin{array}{l}
2^{r_1-r_2}=1\\
 e^{(r_1-r_2)\ln 2}=1= e^{2\pi in}\\
(r_1-r_2)\ln 2=2\pi in\\
r_1-r_2=\sqrt{1-4\la}=2\pi ni/\ln (2)\end{array}
\end{eqnarray}
$$(ref5)
Thus to obtain nontrivial solutions we require $1-4\la <0$ which
implies \framebox{$\la>\frac{1}{4}$}. Thus for $\la >\frac{1}{4}$

$$
\begin{eqnarray}
\sqrt{1-4\la}=i\sqrt{4\la -1}=2\pi ni/\ln (2).
\end{eqnarray}
$$(ref6)
The Eigenvalues are:

$$
\begin{eqnarray}
\la_n =\frac{1}{4}+\frac{\pi^2 n^2}{[\ln (2)]^2},\, 4\la_n
-1=\frac{4\pi^2n^2}{[\ln (2)]^2}=(2\beta_n)^2\quad \beta_n=(n\pi
/\ln 2).
\end{eqnarray}
$$(ref7)
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
$$(ref8)
 Now let us consider expanding a function $f(x)$ in terms of a
`Fourier Series' of these new eigenfunctions in the following form

$$
\begin{eqnarray}
f(x)=\sum\limits_{n=1}^{\infty }c_{n}\phi
_{n}(x)=\sum\limits_{n=1}^{\infty }c_{n}x^{-1/2}\sin \left(
\frac{n\pi \ln x}{\ln 2}\right) \label{eq_genFS}
\end{eqnarray}
$$(ref9)
 In order to determine the coefficients $c_n$ we project the
function $f(x)$ onto the basis functions $\phi_n(x)$ as follows:
\[
\int\limits_{1}^{2}f(x)x^{-1/2}\sin \left( \frac{m\pi \ln x}{\ln
2}\right) dx=\sum\limits_{n=1}^{\infty
}c_{n}\int\limits_{1}^{2}x^{-1/2}\sin \left(
\frac{m\pi \ln x}{\ln 2}\right) x^{-1/2}\sin \left( \frac{n\pi \ln x}{\ln 2}
\right) dx
\]
 The integrals under the summation can be evaluated by making
the simple substitution $z=\ln x$ \ so that $dz=\frac{dx}{x}$, in
which case:
\[
\int\limits_{1}^{2}\sin \left( \frac{m\pi \ln x}{\ln 2}\right) \sin
\left( \frac{n\pi \ln x}{\ln 2}\right)
\frac{dx}{x}=\int\limits_{0}^{\ln 2}\sin \left( \frac{m\pi z}{\ln
2}\right) \sin \left( \frac{n\pi z}{\ln 2}\right) dz=\delta
_{mn}\frac{\ln 2}{2}
\]
 Substituting this result into (\ref{eq_genFS}) we obtain the
following expression for the Fourier coefficients $c_n$:
\[
c_{n}=\frac{2}{\ln 2}\int\limits_{1}^{2}f(x)x^{-1/2}\sin \left(
\frac{n\pi \ln x}{\ln 2}\right) dx
\]

\vfill\eject
\begin{example} A variable coefficient Heat Conduction Problem with Cauchy-Euler Eigenfunctions:

$$
\begin{eqnarray}\begin{array}{ll}
u_t=D{(x^2u_x)}_x-u\quad &1<x<2\quad t>0\\
u(1,t)=0=u(2;t)&u(x,0)=f(x).\end{array}
\end{eqnarray}
$$(ref10)
\end{example}
Let

$$
\begin{eqnarray}
u(x,t)&=&X(x)T(t)\\
\frac{\dot{T}(t)}{DT(t)}&=&\frac{{(x^2X^\prime )}^\prime}{X}-\frac{1}{D}\\
\frac{\dot{T}(t)}{DT(t)}+\frac{1}{D}&=&\frac{{(x^2X^\prime )}^\prime}{X}=-\la
\end{eqnarray}
$$(ref11)

$$
\begin{eqnarray}\begin{array}{ll}
\dot{T}+(1+D\la )T=0 \quad &T(t)=c e^{-(1+D\la )t}\\
\left.\begin{array}{c}{(x^2X^\prime )}^\prime +\la X=0\\
X(1)=0=X(2)\end{array}\right\} \quad &\begin{array}{l}\la_n =\frac{1}{4}+\frac{(\pi n)^2}{[\ln (2)]^2}; X_n(x)=x^{-\frac{1}{2}}\sin\left( n\pi\frac{\ln x}{\ln 2}\right)\\
n=1,2,\ldots \end{array}\end{array}\end{eqnarray}
$$(ref12)

$$
\begin{eqnarray}
u(x,t)&=& x^{-\frac{1}{2}}\sum\limits_{n=1}^\infty c_n e^{-(1+D\la_n)t}\sin\left( n\pi\frac{\ln x}{\ln 2}\right)\\
f(x)&=& u(x,0)=x^{-\frac{1}{2}}\sum\limits_{n=1}^\infty
c_n\sin\left( n\pi\frac{\ln x}{\ln 2}\right)
\end{eqnarray}
$$(ref13)
 where
\[
c_{n}=\frac{2}{\ln 2}\int\limits_{1}^{2}f(x)x^{-1/2}\sin \left(
\frac{n\pi \ln x}{\ln 2}\right) dx
\]

### Solving Laplace'S Equation Using Cauchy-Euler Eigenfunctions
\begin{figure}[htbp]
\begin{center}
\includegraphics[angle=0,width = 6.5cm,clip]{Laplace_Radial_2pt_BVP_Cauchy_Euler_Crop.eps}
\caption{Annular sector subtending an arc of $\alpha$ radians
between the radii $a$ and $b$: $1=b<a=2$} \label{fig:smalldelta}
\end{center}
\end{figure}


$$
\begin{eqnarray}
\Delta u&=& u_{rr}+\frac{1}{r}u_r+\frac{1}{r^2}u_{\theta\theta}=0\quad 1<r<2,\  0<\theta<\al\phantom{\al\al}\\
u(r,0)&=&0,\ u(r,\al )=f( r )\\
u(1,\theta )&=&0, u(2,\theta )=0\\
u(r,\theta )&=&R( r )\Theta (\theta )\\
\frac{r^2R^{\prime\prime}+rR^\prime}{R( r )}&=&
-\frac{\Theta^{\prime\prime}}{\Theta }=-\la^2\quad \mbox{(because of
Homog. BC)}
\end{eqnarray}
$$(ref14)

\noindent $\mathbf{\Theta }$: $\dst\begin{array}{ll}
\Theta^{\prime\prime}-\la^2\Theta =0\quad &\Theta =C\cosh\la\theta +D\sinh\la\theta\\
\Theta (0)=0 &\Theta (0)=C=0\Rightarrow\Theta (\theta
)=D\sinh\la\theta\end{array}$

\noindent $\mathbf{R}$: $r^2R^{\prime\prime}+rR^\prime +\la^2 R=0
(\star), \quad R(1)=0=R(2)$ Although we can easily see that dividing
through by $r$ we can reduce $(\star )$ to S-L form, let us use the
integrating factor

$$
\begin{eqnarray}
\mu(r) =\frac{1}{P} e^{\int\frac{Q}{P}dr}=\frac{1}{r^2}
e^{\int\frac{r}{r^2}dr}=\frac{1}{r^2} e^{\ln r}=\frac{1}{r}.
\end{eqnarray}
$$(ref15)
Therefore $-\frac{1}{r}\cdot (\star )$ and some rearrangement
implies:
\[
-rR^{\prime \prime }-R^{\prime }=-\left( rR^{\prime }\right) ^{\prime }=
\frac{\lambda ^{2}}{r}R
\]

Now let us look for Eigenvalues and Eigenfunctions to $(\star )$.
Let

$$
\begin{eqnarray}
R( r ) &=& r^\gamma\Rightarrow\gamma (\gamma -1)+\gamma +\la^2 =\gamma^2 +\la ^2=0\quad\gamma =\pm i\la .\\
R( r ) &=& c_1r^{i\la }+c_2r^{-i\la}\quad r^{i\la }= e^{i\la \ln r}\\
  &=& A\cos (\la\ln r)+B\sin (\la \ln r)\\
R(1) &=& A\cos\left[\la (\ln 1)\right] +B\sin (\la\ln 1)=A=0\\
R(2) &=& B\sin \left[\la\ln 2\right] =0\Rightarrow\la_n \ln
2=n\pi\quad n=1,2,\ldots
\end{eqnarray}
$$(ref16)
and the corresponding Eigenfunctions are $R_n=\sin\left(
n\pi\frac{\ln r}{\ln 2}\right) $. Therefore

$$
\begin{eqnarray}
u(r,\theta ) =\sum\limits_{n=1}^\infty B_n\sinh(\la _n\theta )\sin
(\la_n \ln r).
\end{eqnarray}
$$(ref17)
Now match BC:

$$
\begin{eqnarray}
f( r )=u(r,\al )=\sum\limits_{n=1}^\infty B_n\sinh (\la_n \al
)\sin\left( n\pi\frac{\ln r}{\ln 2}\right)\,
\end{eqnarray}
$$(ref18)
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
$$(ref19)
Therefore

$$
\begin{eqnarray}
B_n=\frac{2}{\ln 2\sinh\left(\frac{n\pi\al}{\ln
2}\right)}\int\limits_1^2 \frac{f( r )}{r}\sin\left(\frac{n\pi \ln
r}{r}\right)\, dr.
\end{eqnarray}
$$(ref20)

