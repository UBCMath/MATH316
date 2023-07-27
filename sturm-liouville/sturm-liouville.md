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

# Sturm-Liouville Boundary Value Problems

In this lecture we abstract the eigenvalue problems that we have
found so useful thus far for solving the PDEs to a general class of
boundary value problems that share a common set of properties. The
so-called {Sturm-Liouville Problems} define a class of
eigenvalue problems, which include many of the previous problems as
special cases. The $S-L$ Problem helps to identify those assumptions
that are needed to define an eigenvalue problems with the properties
that we require.

```{admonition} Key Concepts
Eigenvalue Problems, Sturm-Liouville
Boundary Value Problems; Robin Boundary conditions.
```

Reference Section: Boyce and Di Prima Section 11.1 and 11.2

## Boundary Value Problems And Sturm-Liouville Theory

### Eigenvalue Problem Summary

- We have seen how useful eigenfunctions are in the solution of
various PDEs.
- The eigenvalue problems we have encountered thus far have been relatively simple.

**I: The Dirichlet Problem:**

$$
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X(0)=0=X(L)
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\lambda _{n}=\frac{n\pi }{L},\ n=1,2,\ldots  \\
X_{n}(x)=\sin \left( \frac{n\pi x}{L}\right)
\end{array}
\right.
$$

**II: The Neumann Problem:**

$$
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X^{\prime }(0)=0=X^{\prime }(L)
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\lambda _{n}=\frac{n\pi }{L},\ n=0,1,2,\ldots  \\
X_{n}(x)=\cos \left( \frac{n\pi x}{L}\right)
\end{array}
\right.
$$

**III: The Periodic Boundary Value Problem:**

$$
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X(-L)=0=X(L) \\
X^{\prime }(-L)=0=X^{\prime }(L)
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\lambda _{n}=\frac{n\pi }{L},\ n=0,1,2,\ldots  \\
X_{n}(x)\in \left\{ 1,\cos \left( \frac{n\pi x}{L}\right) ,\sin \left( \frac{
n\pi x}{L}\right) \right\}
\end{array}
\right.
$$

**IV: Mixed Boundary Value Problem A:**

$$
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X(0)=0=X^{\prime }(L)
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\lambda _{k}=\frac{(2k+1)\pi }{2L},\ k=0,1,2,\ldots  \\
X_{n}(x)=\sin \left( \frac{(2k+1)\pi }{2L}x\right)
\end{array}
\right.
$$

**V: Mixed Boundary Value Problem B:**

$$
\left.
\begin{array}{c}
X^{\prime \prime }+\lambda ^{2}X=0 \\
X^{\prime }(0)=0=X(L)
\end{array}
\right\} \Longrightarrow \left\{
\begin{array}{c}
\lambda _{k}=\frac{(2k+1)\pi }{2L},\ k=0,1,2,\ldots  \\
X_{n}(x)=\cos \left( \frac{(2k+1)\pi }{2L}x\right)
\end{array}
\right.
$$

### The Regular Sturm-Liouville Problem

Consider the the following two-point boundary value problem

$$
\begin{eqnarray}\begin{array}{l}
{\big( p(x)y^\prime \big)}^\prime - q(x)y+\la r(x)y=0\quad 0<x<\ell\\
\al_1 y(0)+\al_2 y^\prime (0)=0\quad \beta_1y(\ell )+\beta_2y^\prime
(\ell )=0\end{array}
\end{eqnarray}
$$(ref0)

where $p$, $p^\prime$, $q$ and $r$ are continuous on $0\leq
x\leq\ell$ and $p(x)\geq 0$ and $r(x)>0$ on $0\leq x\leq\ell$.

We define the Sturm-Liouville eigenvalue problem as:

$$
\begin{eqnarray}
\left.\begin{array}{l}\mathcal{L}y=\la ry\quad\mbox{where}\quad \mathcal{L}y=-{(py^\prime )}^\prime +qy\\
\al_1 y(0)+\al_2y^\prime (0)=0\quad\mbox{and}\quad\beta_1 y(\ell )+\beta_2 y^\prime (\ell )=0\\
p(x)>0\mbox{ and }r(x)>0.\end{array}\right\}\quad \mbox{SL}
\end{eqnarray}
$$(ref1)

````{prf:remark}

1. If $p=1$, $q=0$, $r=1$, $\al_1=1$, $\al_2=0$, $\beta_1=1$, $\beta_2=0$
we obtain Problem~(I) above whereas if $p=1$, $q=0$, $r=1$,
$\al_1=0$, $\al_2=1$, $\beta_1=0$, $\beta_2=1$, we obtain
Problem~(II) above. Notice that the boundary conditions for these
two problems are specified at separate points and are called {\it
{separated BC}}. The periodic BC $X(0)=X(2\pi )$ are not separated
so that Problem~(III) is not technically a SL Problem.

2. If $p>0$ and $r>0$ and $\ell <\infty$ then the SL Problem is said
to be regular. If $p(x)$ or $r(x)$ is zero for some $x$ or the
domain is $[0,\infty )$ then the problem is singular.

3. There is no loss of generality in the so-called self-adjoint form $\mathcal{L}y=-(py^\prime)^\prime+qy$
since it is possible to convert a general 2nd order eigenvalue
problem

  $$
  \begin{eqnarray}
  -P(x)y^{\prime\prime}-Q(x)y^\prime +R(x)y=\la y
  \end{eqnarray}
  $$(ref2)

  to self-adjoint form by multiplying by a suitable integrating factor $\mu (x)$

  $$
  \begin{eqnarray}
  -\mu (x)P(x)y^{\prime\prime}-\mu Q(x)y^\prime +\mu (x)R(x)y=\la\mu
  (x)y \label{eq:SLintegratingFactor}
  \end{eqnarray}
  $$(ref3)
  but expanding the differential operator we obtain

  $$
  \begin{eqnarray}
  \mathcal{L}y=-py^{\prime\prime}-p^\prime y^\prime +qy=\la ry.
  \label{eq:SLexpanded}
  \end{eqnarray}
  $$(ref4)
  Thus comparing (\ref{eq:SLexpanded}) and
  (\ref{eq:SLintegratingFactor}) we can make the following
  identifications: $p=\mu P$ and $p^\prime =\mu Q\Rightarrow p^\prime
  =\mu^\prime P+\mu P^\prime =\mu Q$ which is a linear 1st order ODE
  for $\mu$ with integrating factor
  $exp(\dst\int\frac{P^\prime}{P}-\dst\frac{Q}{P}\, dx)$

  $$
  \begin{eqnarray}
  \mu^\prime +\left(\frac{P^\prime}{P}-\frac{Q}{P}\right)\mu
  =0\Rightarrow
    {\left[ P e^{-\int\frac{Q}{P}\, dx}\mu\right]}^\prime =0\quad  \Rightarrow \framebox{
  $\mu = \dst\frac{e^{\int\frac{Q}{P}\, dx}}{P}$}.
  \end{eqnarray}
  $$(ref5)
````

````{prf:example} Reducing a Boundary Value Problem to SL Form

$$
\begin{eqnarray}\phi^{\prime\prime} +x\phi^\prime +\la\phi =0 \label{eqnonSLfrom}\\
\phi (0)=0=\phi(1) \end{eqnarray}
$$(ref6)

We bring (\ref{eqnonSLfrom}) into SL form by multiplying by the
integrating factor

$$
\begin{eqnarray}\begin{array}{l}
\dst\mu =\frac{1}{P}e^{\int\frac{Q}{P}\, dx}=e^{\int x\, dx}=e^{x^2/2},\quad P(x)=1,\quad Q(x)=x,\quad R(x)=1.\\
\dst e^{x^2/2}\phi^{\prime\prime} +e^{x^2/2}x\phi^\prime +\la e^{x^2/2}\phi =0\\
\dst\hspace{.5in} -{\left(e^{x^2/2}\phi^\prime \right)}^\prime =\la e^{x^2/2}\phi\\
\dst p(x)=e^{x^2/2}\quad r(x)=e^{x^2/2}\end{array}\end{eqnarray}
$$(ref7)
````

````{prf:example} Convert the Equation $-y^{\prime\prime}+x^4y^{\prime}=\la y$ to SL Form

$$
\begin{eqnarray}
P=1,\quad Q=-x^4,\quad \mu = e^{-\int x^4\, dx} &=& e^{-x^5/5}\\
\mbox{Therefore } -e^{-x^5/5}y^{\prime\prime}+e^{-x^5/5}x^4y^\prime &=&\la e^{-x^5/5}\\
-{\big(e^{-x^5/5}y^\prime\big)}^\prime &=& \la e^{-x^5/5}y.
\end{eqnarray}
$$(ref8)
````

## Properties of SL Problems

1. **Eigenvalues**:

   - (a) The eigenvalues $\la$ are all real.
   - (b) There are an $\infty$ \# of eigenvalues $\la_j$ with $\la_1<\la_2 <\ldots <\la_j\ra\infty$ as $j\ra\infty$.
   - (c) $\la_j >0$ provided $\dst\frac{\al_1}{\al_2}<0$, $\dst\frac{\beta_1}{\beta_2}>0$ $q(x)>0$.

2. **Eigenfunctions**: For each $\la_j$ there is an eigenfunction $\phi_j(x)$ that is unique up to a multiplicative const{.} and which
satisfy:

   - (a) $\phi_j(x)$ are real and can be normalized so that $\dst\int\limits_0^\ell r(x)\phi_j^2(x)\, dx=1$.

   - (b) The eigenfunctions corresponding to different eigenvalues are orthogonal with respect to the weight function $r(x)$:

    $$
    \begin{eqnarray}
    \int\limits_0^\ell r(x)\phi_j(x)\phi_k(x)\, dx=0\quad j\ne k.
    \end{eqnarray}
    $$(ref9)

   - (c) $\phi_j(x)$ has exactly $j-1$ zeros on $(0,\ell )$.

3. **Expansion Property**: $\{\phi_j(x)\}$ are complete if $f(x)$ is piecewise smooth then

$$
\begin{eqnarray}\begin{array}{lrcl}
&f(x)&=&\sum\limits_{n=1}^\infty c_n\phi_n(x)\\
\mbox{where }&c_n&=&\frac{\int\limits_0^\ell r(x)f(x)\phi_n(x)\, dx}
{\int\limits_0^\ell r(x)\phi_n^2(x)\, dx}\end{array}\end{eqnarray}
$$(ref10)

````{prf:example} Robin Boundary Conditions
$$
\begin{eqnarray}\begin{array}{ll}
X^{\prime\prime}+\la X=0, &\la =\mu^2\\
X^\prime (0)=h_1X(0), &X^\prime (\ell )=-h_2X(\ell )\end{array}
\end{eqnarray}
$$(ref11)

where $h_1\ge 0$ and $h_2 \ge 0$.

$$
\begin{eqnarray}
X(x) & = & A\cos\mu x+B\sin\mu x\\
X^\prime (x) & = & -A\mu\sin\mu x+B\mu\cos\mu x
\end{eqnarray}
$$(ref12)

**BC 1:** $X^\prime (0)=B\mu =h_1X(0)=h_1 A\quad \Rightarrow A=B\mu /h_1$.
**BC 2:** $X^\prime (\ell )=-A\mu\sin (\mu\ell)+B\mu\cos (\mu\ell )  = -h_2X(\ell )=-h_2[ A\cos\mu\ell +B\sin\mu\ell ]$

$$
\begin{eqnarray} \Rightarrow B\left[ -\frac{\mu^2}{h_1}\sin (\mu\ell )+\mu\cos
(\mu\ell )\right] & = &
   -Bh_2\left[\frac{\mu}{h_1}\cos\mu\ell +\sin\mu\ell\right]
\end{eqnarray}
$$(ref13)

$$
\begin{eqnarray}
B\left\{\left( -\frac{\mu^2}{h_1}+h_2\right)\sin\mu\ell +\left(\mu
+\frac{h_2}{h_1}\mu\right)\cos\mu\ell\right\} =0.
\end{eqnarray}
$$(ref14)

Therefore

$$
\begin{eqnarray}
\tan (\mu\ell )=\left[\frac{\mu (h_1+h_2)}{\mu^2-h_1h_2}\right] .
\end{eqnarray}
$$(ref15)

**Case I:**  $\mathbf{h_1}$ **and** $\mathbf{h_2\ne 0}$

\includegraphics{h_1_eq_1and_h_2_eq_1.eps}

$$
\left.
\begin{array}{c}
X_n=\frac{\mu_n}{h_1} \cos\mu_n x+\sin\mu_n x,  \ \mbox{and}\  \mu_n \sim n \pi/\ell \ \mbox{as}\  n\rightarrow\infty\\
\end{array}
\right.
$$

**Case II:**  $\mathbf{h_1 \ne 0} **and** $\mathbf{h_2 = 0}$

\includegraphics{h_1_eq_1and_h_2_eq_0.eps}

$$
\begin{eqnarray}
X_n & = & \frac{\mu_n}{h_1}\cos\mu_n x+\sin\mu_n x\\
& = & \frac{\cos\mu_n (\ell - x)}{\sin\mu_n \ell}
\end{eqnarray}
$$(ref16)

**Case III:**  $\mathbf{h_1\ra\infty,\quad h_2\ne 0}$

\includegraphics{h_1_eq_100and_h_2_eq_1.eps}

$$
\begin{eqnarray}
X_n & = & \sin (\mu_n x)\\
\mu_n & \sim
&\left[\left(\frac{2n+1}{2}\right)\frac{\pi}{\ell}\right]\quad
n=0,1,2,\ldots \ \mbox{as}\  n\rightarrow\infty
\end{eqnarray}
$$(ref17)
````

## Appendix: Some Proofs For Sturm-Liouville Theory

### Lagrange's Identity

$$
\dst\int\limits_0^\ell
(v\mathcal{L} u-u\mathcal{L} v)\, dx= \left. - p(x)u^\prime
v\right|_0^\ell + \left. p(x)uv^\prime\right|_0^\ell
$$

_Proof_: Let $u$ and $v$ be any sufficiently differentiable
functions, then

$$
\begin{eqnarray}
\int\limits_0^\ell v\mathcal{L} u\, dx &=& \int\limits_0^\ell v\left\{ -{(pu^\prime )}^\prime +qu\right\}\, dx\\
&=& \left. - vpu^\prime\right|_0^\ell +\int\limits_0^\ell u^\prime pv^\prime\, dx+\int\limits_0^\ell uqv\, dx\\
&=& \left. - vpu^\prime\right|_0^\ell + \left. upv^\prime\right|_0^\ell +\int\limits_0^\ell u\left\{ - {(pv^\prime )}^\prime +qv\right\}dx\phantom{dxdx}\\
\mbox{Therefore }\int\limits_0^\ell v\mathcal{L} u\, dx &=&  \left.
 - pvu^\prime\right|_0^\ell +\left. puv^\prime\right|_0^\ell
+\int\limits_0^\ell u\mathcal{L} v\, dx. \hspace{.5in}\square
\end{eqnarray}
$$(ref18)

Now suppose that $u$ and $v$ both satisfy the SL boundary
conditions. I.E.

$$
\begin{eqnarray}
\begin{array}{lcl}\al_1 u(0)+\al_2 u^\prime (0)&=&0\\
\al_1 v(0)+\al_2 v^\prime (0)&=&0\end{array}\quad\begin{array}{lcl}
\beta_1 u(\ell )+\beta_2 u^\prime (\ell ) &=&0\\
\beta_1 v(\ell )+\beta_2 v^\prime (\ell )&=&0\end{array}
\end{eqnarray}
$$(ref19)

then

$$
\begin{eqnarray}
\int\limits_0^\ell v\mathcal{L} u\, dx-\int\limits_0^\ell
u\mathcal{L} v\, dx
  &=& - p(\ell )u^\prime (\ell )v(\ell )+p(\ell )u(\ell )v^\prime (\ell )\\
& &\quad +p(0) u^\prime (0)v(0) -p(0)u(0)v^\prime (0)\\
&=& p(\ell )\left\{ +\frac{\beta_1}{\beta_2}u(\ell )v(\ell )+u(\ell )\left( -\frac{\beta_1}{\beta_2}v(\ell )\right)\right\}\\
& &\quad + p(0)\left\{ -\frac{\al_1}{\al_2} u(0)v(0)-u(0)\left( -\frac{\al_1}{\al_2}v(0)\right)\right\}\phantom{vovovo}\\
&=&0.
\end{eqnarray}
$$(ref20)

Thus $\dst\int\limits_0^\ell v\mathcal{L} u\, dx=\int\limits_0^\ell
u\mathcal{L}v\, dx$ whenever $u$ and $v$ satisfy the SL boundary
condition.

````{prf:observation}
- If $\mathcal{L}$ and BC are such that $\dst\int\limits_0^\ell v\mathcal{L}u\, dx=\int\limits_0^\ell u\mathcal{L}v\, dx$
then $\mathcal{L}$ is said to be **self-adjoint**. 
- Notation: if we define $\dst (f,g)=\int\limits_0^\ell f(x)g(x)\, dx$
then we may write $(v,\mathcal{L}u)=(u,\mathcal{L}v)$.
````

### Proofs using Lagrange's Identity

**(1a) $\la_j$ are real-valued.**

Let $\mathcal{L} y=\la ry$ (1)
$\al_1y(0)+\al_2 y^\prime (0)=0$ $\beta_1y(\ell )+\beta_2 y^\prime
(\ell )=0$. Take the conjugate of (1) $\mathcal{L}\bar{y}=\bar{\la}
r\bar{y}$. By Lagrange's Identity:

$$
\begin{eqnarray}
0&=&(\bar{y},\mathcal{L}y)-(y,\mathcal{L}\bar{y})\\
 &=&(\bar{y},r\la y)-(y,r\bar{\la}\bar{y})\\
 &=&\int\limits_0^\ell\bar{y}(x)r\la y(x)\, dx-\int\limits_0^\ell y(x)r(x)\bar{\la }\bar{y}(x)\, dx\\
 &=&(\la -\bar{\la})\int\limits_0^\ell r(x) {\big| y(x)\big|}^2\, dx
\end{eqnarray}
$$(ref21)

Since $r(x){|y(x)|}^2\geq 0$ it follows that $\la
=\bar{\la}\Rightarrow\la$ is real.

**(1c) $\la_j>0$ provided $\frac{\al_1}{\al_2}<0$, $\frac{\beta_1}{\beta_2}>0$, and $q(x)>0$.**

Consider $\mathcal{L}y=-{(py^\prime )}^\prime +qy=\la ry\,
\mbox{(SL)}$ and multiply (SL) by $y$ and integrate from $0$ to
$\ell$:

$$
\begin{eqnarray}
(y,\mathcal{L}y) &=&\int\limits_0^\ell -{(py^\prime )}^\prime
y+qy^2\, dx
    =\la\int\limits_0^\ell r(x)\big[ y(x)\big]^2\, dx\\
\mbox{Therefore } \la &=& \frac{\int\limits_0^\ell -{(py^\prime
)}^\prime y+qy^2\, dx}{\int\limits_0^\ell ry^2\, dx}
   \quad\mbox{this is known as Rayleigh's Quotient.}\nonumber\\
&=& \frac{[ -py^\prime y]_0^\ell +\int\limits_0^\ell p(y^\prime
)^2+qy^2\, dx}{\int\limits_0^\ell ry^2\, dx} \\
&=& \frac{+p(\ell )\frac{\beta_1}{\beta_2}\big[ y(\ell
)\big]^2-p(0)\frac{\al_1}{\al_2}\big[ y(0)\big]^2+\int\limits_0^\ell
p{(y^\prime )}^2+qy^2\, dx} {\int\limits_0^\ell ry^2\, dx}.
\end{eqnarray}
$$(ref22)

Therefore $\la >0$ since the RHS is all positive.

```{note}
If $q(x)\equiv 0$ and $\al_1=0=\beta_1$ then
with $y^\prime (0)=0=y^\prime (\ell )$ we have nontrivial
eigenfunction $y(x)=1$ and eigenvalue $\la =0$. \vfill\eject
```

**(2b) Eigenfunctions corresponding to different eigenvalues are orthogonal.**

Consider two distinct eigenvalues
$\la_j\ne\la_k$ $\la_j:\mathcal{L}\phi_j=r\la_j\phi_j$ and
$\la_k:\mathcal{L}\phi_k=r\la_k\phi_k$. Then

$$
\begin{eqnarray}
0&=&(\phi_k,\mathcal{L}\phi_j)-(\phi_j,\mathcal{L}\phi_k)\quad\mbox{by Lagrange's Identity}\\
&=& (\phi_k,r\la_j\phi_j)-(\phi_j,r\la_k\phi_k)\\
&=&(\la_j-\la_k)\int\limits_0^\ell r(x)\phi_k(x)\phi_j(x)\, dx
\end{eqnarray}
$$(ref23)
now $\la_j\ne\la_{k}$ implies that

$$
\begin{eqnarray}
\int\limits_0^\ell r(x)\phi_k(x)\phi_j(x)\, dx=0.
\end{eqnarray}
$$(ref24)

**(3) The eigenfunctions form a complete set.**

It is difficult to prove the convergence of
the eigenfunction series expansion for $f(x)$ that is piecewise
smooth. However, if we assume the expansion converges then it is a
simple matter to use orthogonality to determine the coefficients in
the expansion: Let $\dst f(x)=\sum\limits_{n=1}^\infty c_n\phi_n
(x)$.

$$
\begin{eqnarray}
\int\limits_0^\ell f(x)\phi_m (x)r(x)\, dx=\sum\limits_{n=1}^\infty
c_n\int\limits_0^\ell r(x)\phi_m(x)\phi_n(x)\, dx
\end{eqnarray}
$$(ref25)

orthogonality implies

$$
\begin{eqnarray}
c_m=\frac{\int\limits_0^\ell r(x)f(x)\phi_m (x)\,
dx}{\int\limits_0^\ell r(x)\big[\phi_m(x)\big]^2\, dx}.
\end{eqnarray}
$$(ref26)
