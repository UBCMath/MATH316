$$
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

# The Wave Equation On Finite Domains - Solution By Separation Of Variables

In this lecture we discuss the solution of the one dimensional wave
equation on a finite domain using the method of saparation of
variables. The process proceeds in much the same was as with the
heat equation. However, in this case the time equation is a second
order ODE which has an indicial equation with complex roots, which
lead to time functions that are sines and cosines rather than the
exponential decay, which was the case with the heat equation.
Depending on the boundary conditions for the spatial ODE we obtain
the same eigenvalue problems as we did for the case of the heat
equation. Each of these eigensolutions are associated with
particular periodic extension, e.g. the Dirichlet BC  give rise to
eigenfunctions that are sines that are associated with the
odd periodic extension of the solution defined on the domain
$(0,L)$. We will demonstrate, using separation of variables, that
the solution of the wave equation on a finite domain is none other
than the D'Alembert solution in which the initial condition
functions are the periodic extensions of the initial conditions that
correspond to the boundary conditions that apply to the particular
problem.

```{admonition} Key Concepts
The one dimensional Wave Equation;
Finite Domains; Separation of Variables; Even and Odd Extensions and
D'Alembert's solution for finite domains.
```

Reference Section: Boyce and Di Prima Section 10.7

## Solution of the 1D Wave Equation on Finite Domains

### Solution By Separation of Variables

````{prf:example}
$$
\begin{eqnarray}
u_{tt} & = & c^2u_{xx}\qquad 0<x<L,\quad t>0\\
\mbox{BC: } u(0,t) & = & 0,\quad u(L,t)=0\\
\mbox{IC: } u(x,0)& = & f(x),\quad u_t(x,0)=g(x)
\end{eqnarray}
$$(ref0)

For a guitar string $\dst c=\sqrt{\frac{T_0}{\rho_0}}$ whereas for an elastic bar $\dst c=\sqrt{\frac{E}{\rho}}$.

__Separate Variables__  $u(x,t)=X(x)T(t)$

$$
\begin{eqnarray}
\frac{\ddot{T}(t)}{c^2T(t)}=\frac{X^{\prime\prime}(x)}{X(x)}=-\la^2
\end{eqnarray}
$$(ref1)

$$
\begin{eqnarray}
\ddot{T}(t)+\la^2 c^2T(t) =0 & \Rightarrow &T(t)=c_1\cos (\la ct)+c_2\sin (\la ct)\\
\left. \begin{array}{l}X^{\prime\prime}+\la^2 X=0\\
X(0)=0=X(L)\end{array}\right\}
   & \Rightarrow &\left. \begin{array}{l}X(x)=A\cos (\la x)+B\sin\la x  \nonumber\\
        X(0)=A=0\; X(L)=B\sin\la L=0\end{array}\right\}  \nonumber\\
& &\hspace{.55in}\begin{array}{lcl}\la _n&=&\dst\frac{n\pi}{L}\quad n=1,2,\ldots\\
   X_n &=&\sin\left(\dst\frac{n\pi x}{L}\right)\end{array}.
   \nonumber
\end{eqnarray}
$$(ref2)

Therefore

$$
\begin{eqnarray}
u(x,t) & = & \sum\limits_{n=1}^\infty A_n\cos\left(\frac{n\pi
ct}{L}\right)
   \sin\left(\frac{n\pi x}{L}\right) +B_n\sin\left(\frac{n\pi ct}{L}\right)
   \sin\left(\frac{n\pi x}{L}\right)\\
u(x,0) & = & \sum\limits_{n=1}^\infty A_n \sin\left(\frac{n\pi
x}{L}\right)
   =f(x)\Rightarrow\framebox{$A_n=\frac{2}{L}\int\limits_0^L f(x)\sin\left(
   \frac{n\pi x}{L}\right)$}\\
u_t(x,t) & = & \sum\limits_{n=1}^\infty - A_n\left(\frac{n\pi
c}{L}\right)
   \sin\left(\frac{n\pi ct}{L}\right)\sin\left(\frac{n\pi x}{L}\right) +B_n
   \left(\frac{n\pi c}{L}\right)\nonumber\\
& &\hspace{2.5in}\cos\left(\frac{n\pi ct}{L}\right) \sin
   \left(\frac{n\pi x}{L}\right)\\
u_t(x,0) & = & \sum\limits_{n=1}^\infty B_n \left(\frac{n\pi c}
{L}\right) \sin
  \left(\frac{n\pi x}{L}\right) =g(x)\Rightarrow \framebox{$B_n
 \left(\frac{n\pi c}{L}\right) = \frac{2}{L}\int\limits_0^L g(x)
  \sin\left(\frac{n\pi x}{L}\right)\, dx$}\, .\nonumber
\\
\end{eqnarray}
$$(ref3)

Therefore

$$
\begin{eqnarray}
u(x,t)=\sum\limits_{n=1}^\infty \left\{ A_n\cos\left(\frac{n\pi
ct}{L}\right)
   +B_n\sin\left(\frac{n\pi ct}{L}\right)\right\} \sin\left(\frac{n\pi x}{L}\right)
   . \label{eq1DweqFSsol}
\end{eqnarray}
$$(eq1DweqFSsol)
````

````{prf:observation}
1. Period and Frequency of Vibration.

$$
\begin{eqnarray}
\cos\left(\frac{n\pi c}{L}(t+T)\right) =\cos\left(\frac{n\pi
ct}{L}\right)
   \mbox{ provided $\dst\frac{n\pi cT}{L}=2\pi$}
\end{eqnarray}
$$(ref5)

thus $\dst T_n=\left(\frac{2L}{c}\right)\frac{1}{n}$ is the period
(seconds per cycle) of mode $n$. $\dst
f_n=\frac{1}{T_n}=n\left(\frac{c}{2L}\right)$ are the natural
frequencies of vibration.

2. Modes of Vibration: Standing Waves of Wavelength $\la _n=\dst\frac{2L}{n}$.

In the following four figures we plot the fist four modes of vibration. The first, known as 
the fundamental mode of vibration, is associated with the lowest frequency $ \dst
f_1=\frac{1}{T_1}=\left(\frac{c}{2L}\right)$. All higher
frequencies, also known as overtones, are integer multiples of
this fundamental frequency. The nodes  in these modal plots
are indicated by solid circles, which represent the points at which
the displacement associated with a given mode is zero.

```{figure} ../img/wave/vibration_mode_1.png
:name: vibration_mode_1
:alt: vibration_mode_1
:align: center

The fundamental mode of vibration with 2 nodes,
$X_1(x)=\sin\left(\frac{\pi x}{L}\right) \nonumber$.
```

```{figure} ../img/wave/vibration_mode_2.png
:name: vibration_mode_2
:alt: vibration_mode_2
:align: center

The second mode of vibration or first overtone with 3 nodes,
$X_2(x)=\sin\left(\frac{2\pi x}{L}\right) \nonumber$.
```

```{figure} ../img/wave/vibration_mode_3.png
:name: vibration_mode_3
:alt: vibration_mode_3
:align: center

The third mode of vibration with 4 nodes,
$X_3(x)=\sin\left(\frac{3\pi x}{L}\right) \nonumber$.
```

```{figure} ../img/wave/vibration_mode_4.png
:name: vibration_mode_4
:alt: vibration_mode_4
:align: center

The fourth mode of vibration with 5 nodes,
$X_4(x)=\sin\left(\frac{4\pi x}{L}\right) \nonumber$.
```
````

### Interpretation of the Fourier Series Solution in Terms of D'Alembert's Solution

Recall the double angle trigonometric identities

$$
\begin{eqnarray}
\sin (A\pm B)&=&\sin A\cos B  \pm \cos A\sin B; \nonumber\\
   \cos (A\pm B)&=&\cos A\cos \mp \sin A \sin B, \label{doubleangle}
\end{eqnarray}
$$(doubleangle)

which we are going to use to interpret the solution
{eq}`eq1DweqFSsol` in terms of D'Alembert's Solution for an
infinite domain. Using {eq}`doubleangle` we obtain

$$
\begin{eqnarray}
\cos\left(\frac{n\pi ct}{L}\right) \sin\left(\frac{n\pi x}{L}\right)
   & = & \frac{1}{2}\left\{\sin\frac{n\pi}{L} (x+ct)+\sin\left(\frac{n\pi}{L}\right)
   (x-ct)\right\}\phantom{\int\int\int}\\
\sin\left(\frac{n\pi x}{L}\right) \sin\left(\frac{n\pi ct}{L}\right)
   & = & \frac{1}{2}\left\{\cos\frac{n\pi }{L} (x-ct)-\cos\frac{n\pi}{L}
   (x+ct)\right\}
\end{eqnarray}
$$(ref11)

Now

$$
\begin{eqnarray}
\sum\limits_{n=1}^\infty A_n\cos\left(\frac{n\pi ct}{L}\right)
   \sin\left(\frac{n\pi x}{L}\right)
   & = & \frac{1}{2}\sum\limits_{n=1}^\infty A_n
   \left[\sin\left(\frac{n\pi}{L}\right) (x+ct)\right.\nonumber\\
   & &\hspace{.2in}\left. +\sin
   \left(\frac{n\pi}{L}\right) (x-ct)\right]\\
& = & \frac{1}{2}\left[ f_o(x+ct)+f_o(x-ct)\right]
\label{displacementterm}
\end{eqnarray}
$$(displacementterm)

where $f_0$ is the odd periodic extension of $f$. Similarly,

$$
\begin{eqnarray}
\sum\limits_{n=1}^\infty B_n\sin\left(\frac{n\pi ct}{L}\right)
   \sin\left(\frac{n\pi x}{L}\right) =\frac{1}{2}\sum\limits_{n=1}^\infty
   B_n\left[\cos\frac{n\pi}{L} (x-ct)-\cos\frac{n\pi}{L} (x+ct)\right]=\frac{1}{2}\left[G(x-ct)-G(x+ct)\right]\phantom{\int} .
\end{eqnarray}
$$(ref13)

where

$$
\begin{eqnarray}
G(x):=\frac{1}{2}\sum\limits_{n=1}^\infty B_n\cos(\frac{n\pi}{L} x)
\quad \mbox{and} \quad  B_n=\frac{b^g_n}{\lambda_n c L}=\int_0^L
g(x) \sin \lambda_n x dx
\end{eqnarray}
$$(ref14)

$$
\begin{eqnarray*}
G(x) &:&=\sum\limits_{n=1}^{\infty }B_{n}\cos (\frac{n\pi }{L}x)\quad
\mbox{and}\quad B_{n}=\frac{b_{n}^{g}}{\lambda _{n}c}=\frac{2}{\lambda _{n}cL
}\int_{0}^{L}g(x)\sin \lambda _{n}xdx \\
&=&\sum\limits_{n=1}^{\infty }\frac{b_{n}^{g}}{\lambda _{n}c}\cos (\frac{
n\pi }{L}x)\quad \mbox{and}\quad
b_{n}^{g}=\frac{2}{L}\int_{0}^{L}g(x)\sin \lambda _{n}xdx
\end{eqnarray*}
$$(ref15)

therefore

$$
G^{\prime }(x)=-\frac{1}{c}\sum\limits_{n=1}^{\infty }b_{n}^{g}\sin (\frac{
n\pi }{L}x)=-\frac{1}{c}g_{o}(x)
$$

Thus

$$
G(x)=-\frac{1}{c}\int_{0}^{x}g_{o}(s)ds+D
$$

$$
\begin{eqnarray}
\sum\limits_{n=1}^{\infty }B_{n}\sin \left( \frac{n\pi ct}{L}\right)
\sin \left( \frac{n\pi x}{L}\right)  &=&\frac{1}{2}\left[
G(x-ct)-G(x+ct)\right]
\\
&=&\frac{1}{2c}\left\{ \left(
-\int\limits_{0}^{x-ct}g_{o}(s)\,ds+D\right)
-\left( -\int\limits_{0}^{x+ct}g_{o}(s)\,ds+D\right) \right\}  \\
&=&\frac{1}{2c}\int\limits_{x-ct}^{x+ct}g_{o}(s)\,ds
\label{velocityterm}
\end{eqnarray}
$$(velocityterm)

Therefore, combining {eq}`displacementterm` and
{eq}`velocityterm` we obtain the following expression for the
solution of the wave equation for a finite domain in the form of
D'ALembert's solution

$$
\begin{eqnarray}
\framebox{$u(x,t)=\dst\frac{1}{2}\left[ f_o(x+ct)+f_o(x-ct)\right]
    +\dst\frac{1}{2c}\int\limits_{x-ct}^{x+ct} g_o(s)\, ds$}
    \label{eqDAlemGuitar}
\end{eqnarray}
$$(eqDAlemGuitar)

where $f_o$ and $g_o$ are the odd periodic extensions of $f$ and $g$
on $[0,L]$ i.e.

$$
\begin{eqnarray}
f_o(x) & = & \left\{\begin{array}{rrll}
f(x) &0<x<L &\mbox{ and }&f_0(x+2L)=f_0(x_0)\\
-f(-x) &-L<x<0 &&\end{array}\right. \\
g_o(x) & = & \left\{\begin{array}{rrll}
g(x) &0<x<L &\mbox{ and }&g_o(x+2L)=g_o(x)\\
-g(-x) &-L<x<0 &&\end{array}\right. .
\end{eqnarray}
$$(ref18)

````{prf:observation}
1. Equation {eq}`eqDAlemGuitar` above shows that the Wave Equation Solution for a
string tied down at its ends is given by D'Alembert's Solution (see
(23.25) in Lecture 23) in which the initial displacement function is
given by the odd periodic extension $f_0$ of the initial
displacement of the string, and the initial velocity function is
given by the odd periodic extension of $g_0$.

2. Information is carried along the characteristic curves $x+ct=$
const $\quad x-ct=$ const.
\vspace{1in}

3. Observe that the time dependence of the solution involves
$\dst\sin\left(\frac{n\pi ct}{L}\right)$ and
$\cos\left(\dst\frac{n\pi ct}{L}\right)$ which do not decay with
time. Thus the solutions to the Wave Equation persist with time,
whereas the solutions to the Heat Equation typically decay
exponentially with time.
````
