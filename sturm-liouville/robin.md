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

# The Heat Equation With Robin Bc

In this lecture we demonstrate the use of the Sturm-Liouville
eigenfunctions in the solution of the heat equation. We first
discuss the expansion of an arbitrary function $f(x)$ in terms of
the eigenfunctions $\{\phi_n(x)\}$ associated with the Robins
boundary conditions. This is a generalization of the Fourier Series
approach and entails establishing the appropriate normalizing
factors for these eigenfunctions. We then uses the new generalized
Fourier Series to determine a solution to the heat equation when
subject to Robins boundary conditions.

```{admonition} Key Concepts
Eigenvalue Problems, Sturm-Liouville
Boundary Value Problems; Robin Boundary conditions.
```

\vskip 0.1 in
 Reference Section: Boyce and Di Prima Section 11.1 and 11.2

## Solving The Heat Equation With Robin Bc

### Expansion In Robin Eigenfunctions
In this subsection we consider a Robin problem in which $\ell=1$,
$\mathbf{h_1\ra\infty,\ \mbox{and}\  h_2=1}$, which is a Case III
problem as considered in lecture 30. In particular:

\begin{minipage}[t]{4.0in}
\[
\left.
\begin{array}{c}
\phi^{\prime \prime }+\mu^{2}\phi=0 \\
\phi(0)=0, \phi^{\prime }(1)=-\phi(1) \\
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\phi_n  =  \sin (\mu_n x),\\
 \tan(\mu_n)=-\mu_n \\
 \mu_n  \sim
\left[\left(\frac{2n+1}{2}\right){\pi}\right]\ \mbox{as}\
n\rightarrow\infty
\end{array}
\right.
\]
\end{minipage}
\hskip.5in
\begin{minipage}[t!]{2.0in}
\includegraphics{h_1eq_1000and_h_2_eq_1.eps}
\end{minipage}

 Assume that we can expand $f(x)$ in terms of
$\phi_n(x)$:

$$
\begin{eqnarray}
f(x)=\sum\limits_{n=1}^\infty c_n\phi_n(x)
\end{eqnarray}
$$(ref0)

$$
\begin{eqnarray}
\int\limits_0^1 f(x)\sin (\mu_nx)\, dx &=&
c_n\int\limits_0^1\big[\phi_n(x)\big]^2\, dx\\
&=& c_n\frac{1}{2}\big[ 1+\cos^2\mu_n\big]
\end{eqnarray}
$$(ref1)

Therefore

$$
\begin{eqnarray}
c_n=\frac{2}{[1+\cos^2\mu_n ]} \int\limits_0^1 f(x)\sin (\mu_nx)\,
dx.
\end{eqnarray}
$$(ref2)
If $f(x)=x$ then

$$
\begin{eqnarray}
\begin{array}{llcl}
&\int\limits_0^1 x\sin (\mu_nx)\, dx &=
   &{\left. -\frac{\cos (\mu_nx)}{\mu_n}-x\right| }_0^1
      +\frac{1}{\mu_n}\int\limits_0^1\cos\mu_n x\, dx\\
&&= &-\frac{\cos (\mu_n)}{\mu_n}+
   {\left. \frac{\sin\mu_nx}{\mu_n^2}\right| }_0^1\\
\mbox{but }-\mu_n\cos\mu_n =\sin\mu_n&&&\\
&&=
&\frac{\sin\mu_n-\mu_n\cos\mu_n}{\mu_n^2}=2\frac{\sin\mu_n}{\mu_n^2}.
\end{array}
\end{eqnarray}
$$(ref3)
Therefore

$$
\begin{eqnarray}
c_n &=& \frac{4\sin\mu_n}{\mu_n^2 [1+\cos^2\mu_n]}\\
f(x) &=& 4\sum\limits_{n=1}^\infty\frac{\sin\mu_n\sin
(\mu_nx)}{\mu_n^2[1+\cos^2\mu_n]}
\end{eqnarray}
$$(ref4)

### Solving The Heat Equation With Robin Bc

\begin{figure}[htbp]
\begin{center}
\subfigure{\includegraphics[angle=0,width = 3.1
in,clip]{Robin_BC_xt.eps} \label{fig:Robin_xt}} \subfigure[Solution
profiles $u(x,t)$ at various times]
{\includegraphics[angle=0,width=3.1 in,clip]{Robin_BC_uxt.eps}
\label{fig:Robin_uxt}} \caption{Left: Initial and boundary
conditions; Right:Solution profiles $u(x,t)$} \label{fig:RobinHeat}
\end{center}
\end{figure}


$$
\begin{eqnarray}
u_t &=&\al^2 u_{xx}\quad 0<x<1\\
u(0,t)&=&1\quad u_x(1,t)+u(1,t)=0\\
u(x,0)&=& f(x).
\end{eqnarray}
$$(ref5)
Look for a steady state solution $v(x)$

$$
\begin{eqnarray}
\left.\begin{array}{l}v^{\prime\prime}(x)=0\\
v(0)=1\quad v^\prime (1)+v(1)=0\end{array}\right\}
\end{eqnarray}
$$(ref6)

$$
\begin{eqnarray}\begin{array}{c}
v=Ax+B\quad v(0)=B=1\quad v^\prime (x)=A\quad v^\prime (1)+v(1)=A+(A+1)=0\\
A= -1/2\end{array}
\end{eqnarray}
$$(ref7)
Therefore

$$
\begin{eqnarray}
v(x)=1-x/2.
\end{eqnarray}
$$(ref8)
Now let $u(x,t)=v(x)+w(x,t)$

$$
\begin{eqnarray*}
u_t=w_t=\al^2 (v^{\prime\prime}\!\!\!\!\!\nearrow +w_{xx}) &\Rightarrow & w_t=\al^2w_{xx}\\
1=u(0,t)=v(0)+w(0,t)=1+w(0,t) &\Rightarrow & w(0,t)=0
\end{eqnarray*}
$$(ref9)

$$
\begin{eqnarray*}\begin{array}{rcllcl}
0=u_x(1,t)+u(1,t)&=&\left\{ v^\prime (1)+\!\!\!\!\!\nearrow v(1)\right\} &+w_x(1,t)+w(1,t)&\Rightarrow &w_x(1,t)+w(1,t)=0\\
f(x)=u(x,0) &=&v(x)+w(x,0)&&\Rightarrow
&w(x,0)=f(x)-v(x).\end{array}
\end{eqnarray*}
$$(ref10)
Let

$$
\begin{eqnarray}
w(x,t) &=& X(x)T(t)\\
\frac{\dot{T}(t)}{\al^2 T(t)} &=&\frac{X^{\prime\prime}}{X}= -\mu^2\\
T(t) &=& c e^{-\al^2\mu^2t}
\end{eqnarray}
$$(ref11)

$$
\begin{eqnarray}
\left.\begin{array}{l}X^{\prime\prime}+\mu^2X=0\\
X(0)=0\quad X^\prime (1)+X(1)=0\end{array}\right\} \quad
\begin{array}{l}\mbox{The $\mu_n$ are solutions of the transcendental
}\\
\mbox{equation: }\tan\mu_n =-\mu_n.\end{array}
\end{eqnarray}
$$(ref12)

$$
\begin{eqnarray}
X_n(x) &=& \sin (\mu_nx)\\
w(x,t) &=& \sum\limits_{n=1}^\infty c_n e^{-\al^2\mu_n^2t}\sin
(\mu_nx)
\end{eqnarray}
$$(ref13)
where

$$
\begin{eqnarray}
f(x)-v(x) &=& w(x,0)=\sum\limits_{n=1}^\infty c_n\sin (\mu_nx)\\
\Rightarrow c_n &=&\frac{2}{[1+\cos^2\mu_n]}\int\limits_0^1 [f(x)-v(x)]\sin (\mu_n x)\, dx\\
u(x,t) &=& 1-\frac{x}{2}+\sum\limits_{n=1}^\infty c_n e^{-\al^2
\mu_n^2t}\sin (\mu_nx).
\end{eqnarray}
$$(ref14)

