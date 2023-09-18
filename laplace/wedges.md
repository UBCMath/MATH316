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

# Wedges With Cut-Outs, Dirichlet And Neumann Problems On Circular Domains

In this lecture we continue with the solution of Laplace's equation
on circular domains. In particular, a wedge-shaped domain with a
circular cut-out, i.e., a pizza slice with a bite taken out of it,
as well as the Dirichlet and Neumann problems on the interior of a
circle. This leads to the well known Poisson Integral Formula
and an interesting application known as Electrical Impedance
Tomography (EIT).

```{admonition} Key Concepts
Laplace's equation; Circular domains;
Pizza Slice-shaped regions; Dirichlet and Neumann Boundary
Conditions, Poisson Integral Formula; Applications to EIT.
```

Reference Section: Boyce and Di Prima Section 10.8

## Wedges With Cut-Outs, Circles, Holes, And Annuli

### Wedges with Cut-Outs

````{prf:example} A Circular Wedge with a Cut-Out
:label: example-laplace-wedges-0 
```{figure} ../img/laplace/bitten_pizza.png
:name: bitten_pizza 
:align: center

Mixed boundary conditions on a wedge shaped domain with a
cut-out {eq}`eqLaplaceDirichletWedgeCutout`.
```

$$
\begin{eqnarray}
u_{rr}+\frac{1}{r}u_r+\frac{1}{r^2}u_{\theta\theta} =0
\end{eqnarray}
$$(ref-laplace-wedges-0)

$$
\begin{eqnarray}\begin{array}{ll}
u_\theta (r,0)=0&u_\theta (r,\alpha )=0\\
u(b,\theta )=0&u(a,\theta )=f(\theta
)\end{array}\label{eqLaplaceDirichletWedgeCutout}
\end{eqnarray}
$$(eqLaplaceDirichletWedgeCutout)

Let $u(r,\theta )=R( r )\Theta (\theta )$.

$$
\begin{eqnarray}
\frac{r^2(R^{\prime\prime}+\frac{1}{r}R)}{R( r
)}=-\frac{\Theta^{\prime\prime}(\theta )}{\Theta (\theta
)}=\lambda^2\Rightarrow\left\{\begin{array}{l}
r^2R^{\prime\prime}+rR^\prime -\lambda^2 R=0\\
\Theta^{\prime\prime}+\lambda^2 \Theta =0\end{array}\right.
\end{eqnarray}
$$(ref-laplace-wedges-2)

$\mathbf{\Theta \mbox{ equation}\rangle}$

$$
\begin{eqnarray}
& &\left.\begin{array}{l}\Theta^{\prime\prime}+\lambda^2\Theta =0\\
\Theta^\prime (0)=0=\Theta^\prime (\alpha )\end{array}\right\}
\begin{array}{lcl}\Theta &=&A\cos\lambda\theta +B\sin\lambda\theta\\
\Theta^\prime (0)&=&B\lambda =0\Rightarrow B\mbox{ or }\lambda =0,\\
\end{array}\\
& &\quad \begin{array}{lcl}\Theta^\prime &= &-A\lambda\sin\lambda\theta +B\lambda\cos\lambda\theta\\
\Theta^\prime (\alpha ) &= &-A\lambda\sin\lambda\alpha =0,\;
\lambda=\frac{n\pi}{\alpha}\; n=0,1,\ldots
\end{array}
\end{eqnarray}
$$(ref-laplace-wedges-3)

$\mathbf{R \mbox{ equation}\rangle}$ $\mathbf{n=0}$:
${(rR_0^\prime )}^\prime =0$ $rR_0^\prime =B_0$ $R_0=A_0+B_0\ln r$.

Note

$$
\begin{eqnarray}
u_0(b,\theta )=R_0(b)\Theta_0 (\theta )=0\Rightarrow R_0(b)=
   A_0+B_0\ln b=0, A_0=-B_0\ln b.\phantom{\int}
\end{eqnarray}
$$(ref-laplace-wedges-4)

Therefore $R_0=B_0\ln (r/b)$. Choose $B_0=1$.

$\mathbf{n\geq 1}$: $r^2R_n^{\prime\prime}+rR_n^\prime
-\lambda^2 R_n=0$ $R( r )=A_nr^{\lambda_n}+B_nr^{-\lambda_n}$

$$
\begin{eqnarray}
R_n(b) & = & A_nb^{\lambda_n}+B_nb^{-\lambda_n}
   =0\Rightarrow B_n=-A_nb^{2\lambda_n}\\
R_n( r ) & = & A_n[r^{\lambda_n}-b^{2\lambda_n}r^{-\lambda_n}]\quad
   \mbox{ Choose $A_n=1$.}\\
u_n(r,\theta ) & = & \left[ r^{\left(\frac{n\pi}{\alpha}\right)}
   -b^{2\left(\frac{n\pi}{\alpha}\right)} r^{-\left(\frac{n\pi}{\alpha}\right)}\right]
   \cos\left(\frac{n\pi\theta}{\alpha}\right)\\
u_0(r,\theta ) & = & \ln\left(\frac{r}{b}\right) \cdot 1
\end{eqnarray}
$$(ref-laplace-wedges-5)

Therefore

$$
\begin{eqnarray}
u(r,\theta ) & = & c_0\ln\left(\frac{r}{b}\right)
+\sum\limits_{n=1}^\infty c_n\left[
r^{\left(\frac{n\pi}{\alpha}\right)}-b^{\left(\frac{2n\pi}{\alpha}\right)}
r^{-\left(\frac{n\pi}{\alpha}\right)}\right]\cos\left(\frac{n\pi\theta}{\alpha}\right)\\
u(a,\theta ) & = & f(\theta )= 2\frac{\left( c_0\ln\left(
\frac{a}{b}\right)\right)}{2}+\sum\limits_{n=1}^\infty c_n\left[
a^{\left(\frac{n\pi}{\alpha}\right)}
-b^{\left(\frac{2n\pi}{\alpha}\right)}
r^{-\left(\frac{n\pi}{\alpha}\right)}\right]\cos\left(\frac{n\pi\theta}{\alpha}\right)\phantom{\int\int}\\
& &\quad\quad =\frac{a_0}{2}+\sum\limits_{n=1}^\infty
a_n\cos\left(\frac{n\pi\theta}{\alpha}\right) .
\end{eqnarray}
$$(ref-laplace-wedges-6)

Therefore

$$
\begin{eqnarray}
2c_0\ln (a/b)=\frac{2}{\alpha}\int\limits_0^\alpha f(\theta )\, d\theta .
\end{eqnarray}
$$(ref-laplace-wedges-7)

$$
\begin{eqnarray}
c_n & = & \frac{2}{
  \alpha\left[ a^{\left(\frac{n\pi}{\alpha}\right)}
  -b^{\left(\frac{2n\pi}{\alpha}\right)}a^{-\frac{n\pi}{\alpha}}\right] }
\int\limits_0^\alpha f(\theta )\cos\left(\frac{n\pi\theta}{\alpha}\right)\, d\theta\\
c_0 & = & \frac{1}{\alpha\ln (a/b)}\int\limits_0^\alpha f(\theta )\,
d\theta .
\end{eqnarray}
$$(ref-laplace-wedges-8)
````

````{prf:observation}
:label: observation-laplace-wedges-0
In the special case $f(\theta )=1$,
$c_0=\displaystyle\frac{1}{\log (a/b)}$, and $c_n=0$ $n\geq 1$. By Fourier
basis function orthogonality $\int\limits_0^\alpha 1\cdot
\cos\left(\displaystyle\frac{n\pi\theta}{\alpha}\right)\, d\theta =0$ so that
the solution reduces to:

$$
\begin{eqnarray}
u(r,\theta ) = \frac{\log (r/b)}{\log (a/b)}
\end{eqnarray}
$$(ref-laplace-wedges-9)

which is purely radial i.e. has no $\theta$ dependence.
````

### Problems With A Complete Circle As The Boundary

````{prf:example} Dirichlet Problem in the Interior of a Circle
:label: example-laplace-wedges-1 
```{figure} ../img/laplace/dirichlet_circle.png
:name: dirichlet_circle
:align: center

Dirichlet boundary condition prescribed on the on boundary
of a circle {eq}`eqLaplaceDirichletfullCircle`.
```

$$
\begin{eqnarray}
u_{rr}+\frac{1}{r}u_r+\frac{1}{r^2}u_{\theta\theta} =0
  \qquad 0<r<a,\quad 0<\theta <2\pi .
\end{eqnarray}
$$(eqLaplaceDirichletfullCircle)

$$
\begin{eqnarray}\begin{array}{l}
\mbox{BC: }u(a,\theta )=f(\theta )\qquad u(r,\theta )<\infty\quad r\rightarrow 0\\
\mbox{Periodicity } u(\theta +2\pi )=u(\theta
)\quad\mbox{periodic.}\end{array}
\end{eqnarray}
$$(ref-laplace-wedges-11)

Let $u(r,\theta )=R( r )\Theta (\theta )$.

$$
\begin{eqnarray}
\frac{r^2\left( R^{\prime\prime}+\frac{1}{r}R^\prime\right)} {R( r
)}=-\frac{\Theta^{\prime\prime}}{\Theta}=+\lambda^2\Rightarrow
\begin{array}{l}r^2R^{\prime\prime}+rR^\prime -\lambda^2R=0\;\mbox{Euler Eq.}\\
\Theta^{\prime\prime}+\lambda^2\Theta =0\end{array}
\end{eqnarray}
$$(ref-laplace-wedges-12)

$\mathbf{\Theta \mbox{ equation}\rangle}$
$\Theta^{\prime\prime}+\lambda^2\Theta =0$ $\Theta =A\cos (\lambda\theta
)+B\sin (\lambda\theta )$

$$
\begin{eqnarray}
\theta (-\pi ) & = & A\cos (\lambda\pi )-B\sin (\lambda\pi )=\Theta (\pi )=A\cos (\lambda\pi )\\
& &\quad +B\sin (\lambda\pi )\Rightarrow 2B\sin (\lambda\pi )=0\\
\Theta^\prime (\theta ) & = & -A\lambda\sin (\lambda\theta )+B\lambda\cos (\lambda\theta )\\
\Theta^\prime (-\pi ) & = & A\lambda\sin (\lambda\pi )+B\lambda\cos (\lambda\pi
)=\Theta^\prime (\pi )=-A\lambda\sin (\lambda\pi )\\
& &\quad +B\lambda\cos (\lambda\pi )\Rightarrow 2A\lambda\sin (\lambda\pi )=0.
\end{eqnarray}
$$(ref-laplace-wedges-13)

Therefore $\lambda_n=\displaystyle\frac{n\pi}{\pi}=n$; $n=0,1,2,\ldots$
$\Theta_n(\theta )=\cos (n\theta )$ and $\sin (n\theta )$.

$\mathbf{R \mbox{ equation}\rangle}$ $\mathbf{n=0}$:
$\la_0=0\Rightarrow rR_0^{\prime\prime}+R_0^\prime ={(rR_0^\prime
)}^\prime =0$ $R_0=B_0\ln r+A_0$.

$\mathbf{n\geq 1}$: $\mathbf{\lambda_n=n}$:

$$
\begin{eqnarray}\begin{array}{l}
r^2R_n^{\prime\prime}+rR_n^\prime - (n)^2R_n=0\quad\mbox{Euler Eq.}\\
R_n=r^\gamma\Rightarrow\gamma (\gamma -1)+\gamma - n^2=0\quad\gamma =\pm\lambda_n =\pm n\\
R_n=c_nr^{-n}+d_nr^n.\end{array}
\end{eqnarray}
$$(ref-laplace-wedges-14)

Now since $u(r,\theta )<\infty$ as $r\rightarrow 0$ we must exclude
solutions that blow up. Thus $B_0=0$ and $c_n=0$. Therefore

$$
\begin{eqnarray}
u(r,\theta ) & = & \frac{a_0}{2}+\sum\limits_{n=1}^\infty\big\{
  a_n\cos (n\theta )+b_n\sin (n\theta )\big\} r^n \label{eqDirichletSoln}\\
u(a,\theta ) & = &f(\theta )=\frac{a_0}{2} +\sum\limits_{n=1}^\infty
\big[
  a_n\cos (n\theta )+b_n\sin (n\theta )\big] a^n\\
a_0 & = & \frac{1}{\pi}\int\limits_{-\pi}^\pi f(\theta )\,
d\theta\quad a_n=
  \frac{a^{-n}}{\pi}\int\limits_{-\pi}^\pi f(\theta )\cos (n\theta )\, d\theta\\
& &\hspace{1in}b_n=\frac{a^{-n}}{\pi}\int\limits_{-\pi}^\pi f(\theta
)\sin (n\theta )\, d\theta .\end{eqnarray}
$$(eqDirichletSoln)
````

#### The Poisson's Integral Formula

Using the geometric series it is possible to sum the solution
{eq}`eqDirichletSoln` explicitly. We proceed as follows:

$$
\begin{eqnarray}
u(r,\theta ) & = & \frac{a_0}{2}+\sum\limits_{n=1}^\infty \big(
a_n\cos (n\theta )+b_n\sin (n\theta )\big) r^n\\
& = & \frac{1}{2\pi}\int\limits_{-\pi}^\pi f(\phi )d\phi
+\frac{1}{\pi}\sum\limits_{n=1}^\infty
\left(\frac{r}{a}\right)^n\int\limits_{-\pi}^\pi f(\phi )\big[\cos
n\theta\cos n\phi +\sin (n\theta )\sin (n\phi )\big]\, d\phi\nonumber\\
\\
& = & \frac{1}{\pi}\left\{\frac{1}{2}\int\limits_{-\pi}^\pi f(\phi
)\, d\phi
+\sum\limits_{n=1}^\infty\left(\frac{r}{a}\right)^n\int\limits_{-\pi}^\pi
f(\phi )\cos n(\theta -\phi )\, d\phi\right\}\\
& = & \frac{1}{\pi}\int\limits_{-\pi}^\pi f(\phi
)\left\{\frac{1}{2}+\sum\limits_{n=1}^\infty
\left(\frac{r}{a}\right)^n \cos n(\theta - \phi )\right\}\, d\phi
\end{eqnarray}
$$(ref-laplace-wedges-16)

Now

$$
\begin{eqnarray}\begin{array}{l}
\sum\limits_{n=1}^\infty \left(\frac{r}{a}\right)^n\cos n(\theta
-\phi )\\
\phantom{\frac{r}{a}}\\
\\
\\
\\
\\
\\\end{array}\begin{array}{cl}
=&Re\sum\limits_{n=1}^\infty z^n\\
=&Re\left(\frac{z}{1-z}\right)\\
=&Re\left(\frac{z(1-\bar{z})}{(1-z)(1-\bar{z})}\right)\\
=&\frac{\left(\frac{r}{a}\right)\cos (\theta -\phi
)-\left(\frac{r}{a}\right)^2}
  {1-2\left(\frac{r}{a}\right)\cos (\theta -\phi )+\left(\frac{r}{a}\right)^2}\\
=&\frac{ar\cos (\theta -\phi )-r^2} {a^2-2ar\cos (\theta
  -\phi )+r^2}\end{array}\!\!\!\!\! \begin{array}{rcl}
z&=&\left(\frac{r}{a}\right)e^{i(\theta -\phi )}\\
(1-z)(1-\bar{z})&=&1-(z+\bar{z})+|z|^2\\
&=&1-2\left(\frac{r}{a}\right) \cos (\theta -\phi
)+\left(\frac{r}{a}\right)^2\\
z(1-\bar{z})&=&\frac{r}{a}e^{i(\theta -\phi
)}-\left(\frac{r}{a}\right)^2\\
\\
\\
\\\end{array}\nonumber\end{eqnarray}
$$(ref-laplace-wedges-17)

Therefore

$$
\begin{eqnarray}
u(r,\theta ) &=& \frac{1}{\pi}\int\limits_{-\pi}^\pi f(\phi )
\left\{\frac{1}{2}+\frac{ar\cos (\theta -\phi )-r^2}
  {a^2-2ar\cos (\theta -\phi )+r^2}\right\}\, d\phi\\
&=& \frac{1}{\pi}\int\limits_{-\pi}^\pi f(\phi )\left\{
  \frac{\frac{1}{2}a^2-ar\cos (\theta -\phi )+\frac{1}{2}r^2+ar\cos (\theta -\phi )-r^2}
  {a^2-2ar\cos (\theta -\phi )+r^2}\right\}\, d\phi\nonumber\\
&&
\end{eqnarray}
$$(ref-laplace-wedges-18)

$$
\begin{eqnarray}
\framebox{$\displaystyle u(r,\theta
)=\frac{1}{2\pi}(a^2-r^2)\int\limits_{-\pi}^\pi \frac{f(\phi
)}{a^2-2ar\cos (\theta -\phi )+r^2}\, d\phi $}
\end{eqnarray}
$$(ref-laplace-wedges-19)

#### Domain Exterior to a Circle

For problem exterior to a circle we require that $u(r,\theta
)<\infty$ as $r\rightarrow\infty$. In this case we require that $B_0=0$ and
that $d_n=0$ so that $R_0=A_0$ and $R_n=r^{-n}\cdot\big[ a_n\cos
n\theta +b_n\sin (n\theta )\big]$. In this case the appropriate
solution to the Dirichlet problem becomes:

$$
\begin{eqnarray}
u(r,\theta ) & = & \frac{a_0}{2}+\sum\limits_{n=1}^\infty\big[ a_n\cos (n\theta )+b_n\sin (n\theta )\big] r^{-n}\\
a_n & = & \frac{a^n}{\pi}\int\limits_{-\pi}^\pi f(\theta )\cos
(n\theta )\, d\theta ,b_n=\frac{a^n}{\pi}\int\limits_{-\pi}^\pi
f(\theta )\sin (n\theta )\, d\theta .
\end{eqnarray}
$$(ref-laplace-wedges-20)

````{prf:example} Neumann Problem on the Interior of a Circle
:label: example-laplace-wedges-2 
```{figure} ../img/laplace/neumann_circle.png
:name: neumann_circle
:align: center

Neumann boundary condition prescribed on the on boundary of
a circle {eq}`eqLaplaceNeumannfullCircle`.
```

$$
\begin{eqnarray}\begin{array}{l}
u_{rr}+\displaystyle\frac{1}{r}u_r+\displaystyle\frac{1}{r^2}u_{\theta\theta}=0\\
\displaystyle\frac{\partial u}{\partial r}(a,\theta )=f(\theta )\\
u\quad 2\pi \mbox{ - periodic.}\end{array}
\label{eqLaplaceNeumannfullCircle}
\end{eqnarray}
$$(eqLaplaceNeumannfullCircle)

$$
\begin{eqnarray}
u(r,\theta )= &\displaystyle\frac{a_0}{2}& +\sum\limits_{n=1}^\infty r^n
  \big[ a_n\cos (n\theta )+b_n\sin (n\theta )\big]\\
{\left.\frac{\partial u}{\partial r}\right|}_{r=a}= f(\theta ) & = &
   \sum\limits_{n=1}^\infty nr^{n-1} {\left[ a_n\cos (n\theta )+
  b_n\sin (n\theta )\big]\right|}_{r=a}\\
& = & \sum\limits_{n=1}^\infty na^{n-1}\big[ a_n\cos (n\theta
)+b_n\sin (n\theta )\big] .
\end{eqnarray}
$$(ref-laplace-wedges-22)

A solution will not exist unless
$a_0=\displaystyle\frac{1}{2\pi}\int\limits_{-\pi}^\pi f(\theta )\, d\theta
=0$. Otherwise there is a net flux of heat across the boundary and
no steady state solution will exist.
````

#### Special Case - Electrical Impedance Tomography (EIT)

Assume that $f(\theta )=I_0\delta \left(\theta -\displaystyle\frac{\pi}{2} \right) -
I_0\delta\left(\theta +\displaystyle\frac{\pi}{2}\right)$.

```{figure} ../img/laplace/EIT_BVP.png
:name: EIT_BVP
:align: center

(left) Level sets of potential field {eq}`EITSolution`.
(middle) BVP for $U_1$ field and voltage measurements $\delta U_1$ 
with $\sigma=const$. Actual voltage $\delta V_1$ with a feature (red).
(right) Field $U_0$ used in the sensitivity Thm.

Sensitivity Thm: $\delta V_{1}-\delta U_{1}\simeq
-\frac{1}{I}\int\limits_{\Omega }\nabla U_{1}\cdot \nabla
U_{0}\delta \sigma dv$ \& backprojection rule: $\delta \sigma \simeq
-\frac{\left( \delta V_{1}-\delta U_{1}\right) }{\delta
U_{1}}\sigma$.
```

$$
\begin{eqnarray}
f(-\theta ) & = & I_0\delta \left( -\left(\theta
+\frac{\pi}{2}\right)\right)
   - I_0\delta\left( -\left(\theta -\frac{\pi}{2}\right)\right)\\
& = & I_0\delta\left(\theta +\frac{\pi}{2}\right) -
I_0\delta\left(\theta -\frac{\pi}{2}\right) =-f(\theta )
\end{eqnarray}
$$(ref-laplace-wedges-23)

Thus $f$ is odd $\Rightarrow a_0=a_n=0$.

$$
\begin{eqnarray}
na^{n-1}b_n & = & \frac{2}{\pi}\int\limits_0^\pi I_0\delta\left(\theta -\frac{\pi}{2}\right)\sin (n\theta )\, d\theta\\
b_n & = & \frac{2I_0}{\pi na^{n-1}}\sin\left(\frac{n\pi}{2}\right)\\
u(r,\theta ) & = &
\frac{2aI_0}{\pi}\sum\limits_{n=1}^\infty\frac{\sin (n\theta
)}{n}\sin\left(\frac{n\pi}{2}\right) {\left(\frac{r}{a}\right)}^n
\end{eqnarray}
$$(ref-laplace-wedges-24)

$$
\begin{eqnarray}
\mbox{For enrichment }\downarrow & = &
\frac{2aI_0}{\pi}\sum\limits_{n=1}^\infty\left(\frac{r}{a}\right)^n\frac{1}{2}
\left[\frac{\cos n\left(\theta -\frac{\pi}{2}\right)}{n}-
  \frac{\cos n\left(\theta +\frac{\pi}{2}\right)}{n}\right]\phantom{\int\int\int}\\
& = &
\frac{aI_0}{\pi}\sum\limits_{n=1}^\infty\left(\frac{r}{a}\right)^n
\left[\frac{\cos n\left(\theta -\frac{\pi}{2}\right)}{n}-
   \frac{\cos n\left(\theta +\frac{\pi}{2}\right)}{n}\right]\phantom{\int\int}\\
& = & \frac{aI_0}{\pi}R\mathrm{e} \left[\sum\limits_{n=1}^\infty
{\frac{z_1^n}{n}}-\sum\limits_{n=1}^\infty\frac{z_2^n}{n}\right]
\begin{array}{lcl}
z_1&=&\left(\frac{r}{a}\right)\mathrm{e}^{i\left(\theta -\frac{\pi}{2}\right)}\\
z_2&=&\left(\frac{r}{a}\right)\mathrm{e}^{i\left(\theta
+\frac{\pi}{2}\right)}
\end{array}.\phantom{\int\int}
\end{eqnarray}
$$(ref-laplace-wedges-25)

Now for $\mathbf{|z|<1}$:

$$
\begin{eqnarray}\begin{array}{rcrcl}
\displaystyle\frac{1}{1-z}&=&1+z+z^2+\cdots&=&\displaystyle\sum\limits_{k=0}^\infty z^k\\
-\ln (1-z)&=&z+\displaystyle\frac{z^2}{2}+\cdots
&=&\displaystyle\sum\limits_{n=1}^\infty \displaystyle\frac{z^n}{n}\end{array}.
\end{eqnarray}
$$(ref-laplace-wedges-26)

Therefore

$$
\begin{eqnarray}
u(r,\theta )=-\frac{aI_0}{\pi}R\mathrm{e} \left[
\ln\left(\frac{1-z_1}{1-z_2}\right)\right] .
\end{eqnarray}
$$(ref-laplace-wedges-27)

Now if $(1-z)=A\mathrm{e}^{i\phi}$ then

$$
\begin{eqnarray}
R\mathrm{e} \big[\ln (1-z)\big] =R\mathrm{e} \left[\ln \left(
A\mathrm{e}^{i\phi}\right)\right] =R\mathrm{e} [\ln A+i\phi ]=\ln A.
\end{eqnarray}
$$(ref-laplace-wedges-28)

Therefore

$$
\begin{eqnarray}
u(r,\theta )=-\frac{aI_0}{2\pi}\ln
{\left|\frac{1-z_1}{1-z_2}\right|}^2.
\end{eqnarray}
$$(ref-laplace-wedges-29)

Now

$$
\begin{eqnarray}
z_1 & = & \left(\frac{r}{a}\right) \mathrm{e}^{i\left(\theta -\frac{\pi}{2}\right)}=\rho \mathrm{e}^{i\phi_1}\\
\|1 -z_1\|^2 & = & (1-z_1)\overline{(1-z_1)}=\big( 1-\rho \mathrm{e}^{i\phi_1}\big)\big( 1-\rho \mathrm{e}^{-i\phi_1}\big)\\
& = & 1-\rho\left( \mathrm{e}^{i\phi_1}+\mathrm{e}^{-i\phi_1}\right) +\rho^2\\
& = & 1-2\rho\cos\phi_1 +\rho^2 .
\end{eqnarray}
$$(ref-laplace-wedges-30)

Similarly $z_2=\left(\displaystyle\frac{r}{a}\right) \mathrm{e}^{i\left(\phi
+\frac{\pi}{2}\right)}=\rho \mathrm{e}^{i\phi_2}$ and
$|1-z_2|^2=1-2\rho\cos\phi_2 +\rho^2$. Therefore

$$
\begin{eqnarray}
u(r,\theta ) & = &
\frac{aI_0}{2\pi}\ln\left[\frac{1-2(\frac{r}{a})\cos (\theta
+\frac{\pi}{2}) +{(\frac{r}{a})}^2}
   {1-2(\frac{r}{a})\cos (\theta -\frac{\pi}{2})+{(\frac{r}{a})}^2}\right]\\
u(r,\theta ) & = & \frac{a
I_0}{2\pi}\ln\left[\frac{a^2+2ar\sin\theta +r^2}
   {a^2-2ar\sin\theta +r^2}\right] \label{eq:EITSolution}
\end{eqnarray}
$$(EITSolution)

#### Electrical Impedance Tomography - from the US Patent

```{figure} ../img/laplace/EIT_patent.png
:name: EIT_patent
:align: center

(left) Level sets of potential field used in the US Patent.
(middle) Medical imaging device.
(right) Actual image obtained from EIT.

Figures taken from a US Patent for a medical imaging device
based in EIT and the level sets of the solution given in
{eq}`EITSolution`.
```
