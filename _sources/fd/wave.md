# Wave Equation

Consider the following initial boundary value problem for the wave equation:

$$
\begin{eqnarray}
u_{tt}&=&c^2u_{xx}\quad 0<x<L\\
\mbox{BC:}\quad u(0,t) &=& 0\quad u(L,t)=0\\
\mbox{IC:}\quad u(x,0)&=&f(x)\\
\frac{\partial u}{\partial t}(x,0)&=&g(x) \label{eqwaveVelIC}
\end{eqnarray}
$$(eqwaveVelIC)

We introduce a finite difference mesh $x_n=n\Delta t$, $t_k=k\Delta t$ and
let the corresponding nodal values be denoted by

$$
\begin{eqnarray*}
u_n^k\simeq u(x_n,t_k).
\end{eqnarray*}
$$(ref-fd-pdes-31)

```{figure} ../img/fd/wave_eq_scheme.png
:name: fd_scheme_wave
:alt: fd_scheme_wave
:align: center
```

Now approximating derivatives by central differences both in space
and time we obtain

$$
\begin{eqnarray}
\frac{u_n^{k+1}-2u_n^k+u_n^{k-1}}{\Delta t^2}&=&c^2\left(\frac{u_{n+1}^k-2u_n^k+u_{n-1}^k}{\Delta x^2}\right) +O(\Delta x^2,\Delta t^2).\\
\mbox{Therefore}\quad u_n^{k+1}&=&2u_n^k-u_n^{k-1}+{\left(\frac{c\Delta t}{\Delta x}\right)}^2\big( u_{n+1}^k-2u_n^k+u_{n-1}^k\big)\nonumber\\
\underbrace{u_n^{k+1}}_{\mbox{time level $k+1$}}&=&
\underbrace{r^2u_{n+1}^k+2(1-r^2)u_n^k+r^2u_{n-1}^k}_{\mbox{time
level $k$}}-\underbrace{u_n^{k-1}}_{\mbox{time level $k-1$}}
\label{eq_FD_Wavekgen}
\end{eqnarray}
$$(eq_FD_Wavekgen)

Here $r=(c\Delta t/\Delta x)$ is known as the Courant Number. We observe
that the Discrete Equation {eq}`eq_FD_Wavekgen` involves three
distinct levels in which known data is transferred from steps $k-1$
and $k$ to step $k+1$.

___Initial Conditions - Starting the Solution___

The 3-level scheme poses some challenges when imposing the initial
conditions. If we imagine a row of false mesh points at time $t=-\Delta
t=t_{-1}$, then the initial velocity condition {eq}`eqwaveVelIC`
can be approximated using central differences as:

$$
\begin{eqnarray}
\frac{u_n^1-u_n^{-1}}{2\Delta t}& = & g(x_n)\label{FalseICWaveEq}
\end{eqnarray}
$$(ref-fd-pdes-33)

therefore

$$
\begin{eqnarray}
u_n^{-1}&=&u_n^1-2 \Delta t g(x_n) \label{eqWaveIC}
\end{eqnarray}
$$(eqWaveIC)

Now we assume that the Discrete Wave Equation {eq}`eq_FD_Wavekgen`
also holds at $t=0$ so that

$$
\begin{eqnarray}
u_n^1=r^2u_{n+1}^0+2(1-r^2)u_n^0+r^2u_{n-1}^0-u_n^{-1}
\label{eq_FD_Wavekeq1}
\end{eqnarray}
$$(eq_FD_Wavekeq1)

Substituting {eq}`eqWaveIC` into {eq}`eq_FD_Wavekeq1` and
re-arranging we obtain:

$$
\begin{eqnarray}
u_n^1=\frac{1}{2}(r^2u_{n+1}^0+2(1-r^2)u_n^0+r^2u_{n-1}^0)+\Delta t
g(x_n) \label{eq_FD_Wavekeq1c}
\end{eqnarray}
$$(eq_FD_Wavekeq1c)

Since $u_n^0=f(x_n)$ and $g(x_n)$ are known, we are now in a
position to specify the first two rows of nodes. This is sufficient
to start the recursion {eq}`eq_FD_Wavekgen` for all subsequent
steps.

```{figure} ../img/fd/wave_scheme_false_mesh.png
:name: false_mesh
:alt: wave_scheme_false_mesh
:align: center

8.20 refers to {eq}`eq_FD_Wavekgen` and 8.24 refers to {eq}`eq_FD_Wavekeq1c`.
$\bullet=$ nodes and $\circ=$ false mesh points used to derive {eq}`eq_FD_Wavekeq1c`
but not actually used in the computation.
```

#### Stability of the Finite Difference Scheme for the Wave Equation

Consider the following finite difference approximation to the 1D
wave equation:

$$
\begin{eqnarray}
u_n^{k+1}=r^2u_{n+1}^k+2(1-r^2)u_n^k+r^2u_{n-1}^k-u_n^{k-1}
\label{FDWaveEqGen}
\end{eqnarray}
$$(FDWaveEqGen)

and, as in the case of the heat equation, substitute $\displaystyle
u_n^k=\phi_ke^{in\Delta x\theta}$ into {eq}`FDWaveEqGen`

$$
\begin{eqnarray*}
e^{in\Delta x\theta }\phi _{k+1}=\left( r^{2}e^{i\Delta x\theta
}+2\left( 1-r^{2}\right) +r^{2}e^{-i\Delta x\theta }\right)
e^{in\Delta x\theta }\phi _{k}-e^{in\Delta x\theta }\phi _{k-1}
\end{eqnarray*}
$$(ref-fd-pdes-38)

Canceling terms and using the double angle formulae

$$
\begin{eqnarray}
\phi _{k+1} &=&2\left( 1+r^{2}\left( \cos \Delta x\theta -1\right)
\right)
\phi _{k}-\phi _{k-1} \\
&=&2\left( 1-2r^{2}\sin ^{2}\frac{\Delta x\theta }{2}\right) \phi
_{k}-\phi _{k-1} \label{waveStabEq}
\end{eqnarray}
$$(waveStabEq)

 If we now assume that $\phi _{k}$ has the following exponential
form $\phi _{k}=G^{k}$ then {eq}`waveStabEq` reduces to the
following quadratic equation:

$$
\begin{eqnarray}
G^{2}-2\gamma G+1=0 \label{StabQuadratic}
\end{eqnarray}
$$(ref-fd-pdes-40)

where $\gamma =\left( 1-2r^{2}\sin ^{2}\frac{\Delta x\theta
}{2}\right) $. The solutions of this quadratic equation are given by

$$
\begin{eqnarray}
G_{1,2}=\gamma \pm \sqrt{\gamma ^{2}-1}
\end{eqnarray}
$$(ref-fd-pdes-41)

Now since $G_1$ and $G_2$ are the roots of this quadratic we may
conclude that

$$
\begin{eqnarray}
(G-G_{1})(G-G_{2})=G^{2}-(G_{1}+G_{2})G+G_{1}G_{2}=0
\label{StabQuadratic1}
\end{eqnarray}
$$(StabQuadratic1)

Comparing the last terms in these two quadratic equations
(\ref{StabQuadratic}) and {eq}`StabQuadratic1` we conclude that

$$
\begin{eqnarray}
G_{1}G_{2}=1. \label{unitDiskConstraint}
\end{eqnarray}
$$(unitDiskConstraint)

However, for stability of solutions for the form $\phi _{k}=G^{k}$,
we require that $|G_1|\le 1$ and $|G_2|\le 1$. Given the constraint
{eq}`unitDiskConstraint`, the only possibility, if the solutions
are to be stable, is that $|G_1|=|G_2|= 1$. Thus $G$ must fall on
the unit disk, which implies that

$$
\begin{eqnarray*}
|\gamma|\le 1
\end{eqnarray*}
$$(ref-fd-pdes-44)

Thus

$$
\begin{eqnarray*}
\left| 1-2r^{2}\sin ^{2}\frac{\Delta x\theta }{2}\right|\le 1
\end{eqnarray*}
$$(ref-fd-pdes-45)

or

$$
\begin{eqnarray*}
-1\le  1-2r^{2}\sin ^{2}\frac{\Delta x\theta }{2}\le 1
\end{eqnarray*}
$$(ref-fd-pdes-46)

so that

$$
\begin{eqnarray}
-2\le -2r^{2}\sin ^{2}\frac{\Delta x\theta }{2}\le 0
\label{StabCondWaveEq}
\end{eqnarray}
$$(StabCondWaveEq)

The second inequality in {eq}`StabCondWaveEq` is satisfied
automatically, while the first leads to the condition

$$
\begin{eqnarray*}
r^{2}\sin ^{2}\frac{\Delta x\theta }{2}\le 1
\end{eqnarray*}
$$(ref-fd-pdes-48)

Since the maximum value that $\sin ^{2}(\frac{\Delta x\theta }{2})$
can achieve is 1, we conclude that the condition for stability is

$$
\begin{eqnarray*}
r=(c\Delta t/\Delta x)\le 1
\end{eqnarray*}
$$(ref-fd-pdes-49)

or

$$
\begin{eqnarray}
\Delta t\le  \frac{\Delta x }{c} \label{CFL_Wave}
\end{eqnarray}
$$(CFL_Wave)

The condition {eq}`CFL_Wave`, which imposes an upper bound on the
time step that can be used, is known as the Courant-Friedrichs-Lewy
or CFL condition.