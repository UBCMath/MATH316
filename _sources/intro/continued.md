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

# Introduction To Partial Differential Equations (Continued)

In this lecture we will continue with the derivation of the basic PDEs studied in this course 
from different physical situations. We deliberately explore the different
paths to arrive at the same partial differential equations to
emphasize way in which the models from disparate applications share
the same features.

```{admonition} Key Concepts
Partial Differential Equations (PDEs);
Elliptic, Parabolic, Hyperbolic PDEs; The heat Equation, The Wave
Equation, and Laplace's Equation, Modeling and Derivation of PDEs.
```

## Introduction To PDEs (Continued)

### The Wave Equation

Consider an elastic rod having a density $\rho$ and cross-sectional
area $A$, and let $\sigma (x,t)$ be the pressure in the rod at $x$
at time $t$ and $u(x,t)$ the displacement of the rod from its
equilibrium position.

```{figure} ../img/intro/balanced_momentum.png
:name: balanced_momentum
:alt: balance_of_linear_momentum
:align: center

We consider the net force $\sigma (x+\Delta x,t)A-\sigma
(x,t)A$ acting on the bar, which according to Newton's Second Law of
motion, must be balanced by the product of mass of the bar and its
acceleration
```

Balance of Linear Momentum $"F=Ma"$.

$$
\begin{eqnarray}
\sigma (x+\Delta x,t)A-\sigma (x,t)A &=& \rho A\Delta
x\frac{\partial^2u}{\partial t^2}\nonumber\\
\frac{\sigma (x+\Delta x,t)-\sigma (x,t)}{\Delta x}&=& \rho
\frac{\partial^2 u}{\partial t^2} \nonumber\\
\mbox{Let} \ \Delta x\rightarrow 0, \mbox{which
yields}\quad\frac{\partial\sigma}{\partial x} &=& \rho \frac{\partial^2 u}{\partial t^2}
\label{BLM}
\end{eqnarray}
$$(BLM)

We observe that {eq}`BLM` involves two unknown
quantities the stress in the bar $\sigma$ and the resulting
displacement $u$. In order to have sufficient information to solve
for the unknowns we need an additional equation, which is provided
by a constitutive relation known as Hooke's Law (see figure
{numref}`hookes_law`). Experimental data characterizes the
"stiffness" of the material by the parameter $E$ known as the
Young's Modulus, which  provides  a linear relationship between the
stress to which the bar is subjected and the relative displacement
$\frac{\Delta u}{\Delta x}=\frac{u(x+\Delta x,t)-u(x,t)}{\Delta x}
\approx \frac{\partial u}{\partial x} := \epsilon $, or strain $\epsilon$.

```{figure} ../img/intro/hookes_law.png
:name: hookes_law
:alt: hookes_law
:align: center

The stress on the bar $\sigma$ is related to the strain
$\epsilon$ by Hooke's Law
```

Substituting the stress strain relationship $\sigma =E\frac{\partial u}{\partial
x}$ into {eq}`BLM` we obtain the second order wave equation.

$$
\begin{equation}
\frac{\partial^2u}{\partial t^2}=\left(\frac{E}{\rho }\right)\frac{\partial^2u}{\partial
x^2}=c^2\frac{\partial^2u}{\partial x^2}, \mbox{where} \quad
c=\sqrt{\frac{E}{\rho}} \label{eq:acoustic_wave_eq}
\end{equation}
$$(acoustic_wave_eq)

````{prf:observation}
:label: observation-intro-continued-0
We note that {eq}`acoustic_wave_eq`
is in precisely the same form as the wave equation that occurred in
the context of modeling traffic flow down a highway. Equation
{eq}`acoustic_wave_eq` describes the motion of an acoustic wave
that travels in an elastic bar.
````

### The Heat/Diffusion Equation

```{figure} ../img/intro/heat_flux.png
:name: heat_flux
:alt: heat_flux
:align: center

Heat flux in a conducting bar that occupies the region
$[x,x+\Delta x]$.
```

#### Fourier's Law and Heat Conduction

How do we build a model of the flow of heat in a conductor? 
Consider an elemental length of a conducting
metal bar (think copper or aluminium) that occupies the interval
$[x,x+\Delta x]$. We define the following material properties of the
conductor:

*Heat Capacity:* Let $C$ be the specific heat capacity of the
material, which is defined to be the amount of heat required to
change one kilogram of the material by one degree Kelvin, i.e.,
$[C]=\frac{J}{kg\dot K}$.

*Density:* Let $\rho$ be the density of the conductor,
i.e., $[\rho]=\frac{kg}{m^3}$.

*Heat Flux:* Let $q$ be the flux of heat energy per unit area,
i.e., $[q]=\frac{J}{m^2 s}$.

If $A$ is the cross sectional area of the bar, then the
change of temperature $u(x,t)$ within the element of length $\Delta
x$ over a time interval $\Delta t$ is given by

$$
\begin{equation}
 C \Delta u \rho \Delta x A = C \{ u(x,t+\Delta t)-u(x,t)\}\Delta x A\simeq\{
q(x,t)-q(x+\Delta x,t)\} A \Delta t  \label{eq:changeinHeat}
\end{equation}
$$(changeinHeat)

Now dividing by $A \Delta x \Delta t$ and taking the limit as
$\Delta x, \Delta t \rightarrow 0$, we obtain

$$
\begin{equation}
 \rho C \frac{\partial u}{\partial t}+\frac{\partial q}{\partial x}=0
\label{eq:HeatConservationLaw}
\end{equation}
$$(ref-intro-continued-3)

```{figure} ../img/intro/fourier_law.png
:name: fourier_law
:alt: fourier_law
:align: center

Fourier's Law of Heat Conduction: heat moves from points at
which the temperature is higher in the direction of points at a
lower temperature and the flux is given by $q =-k \frac{\partial u}{\partial x}$
```

It is our common experience that heat flows from hotter regions to cooler ones, which
is illustrated in figure {numref}`fourier_law` in which the
directions of the flux of heat $q$ are indicated depending upon the
sign of $ \frac{\partial u}{\partial x}$. In addition, experimental evidence
suggests that the flux of heat is proportional to the negative of
the spatial gradient of the temperature, which relationship is known
as Fourier's Law of heat conduction

$$
\begin{equation}
 q =-k \frac{\partial u}{\partial x}
\label{eq:FouriersLaw}
\end{equation}
$$(FouriersLaw)

where $k$ is the thermal conductivity having dimensions $[k]=
\frac{J}{s m K}$. Substituting {eq}`FouriersLaw` into
{eq}`changeinHeat` and dividing by $\rho C$ we obtain the heat
equation

$$
\begin{equation}
\frac{\partial u}{\partial t}=\alpha^2\frac{\partial^2u}{\partial x^2}
\label{eq:Heat_equation}
\end{equation}
$$(Heat_equation)

where $ \alpha^2 = \frac{k}{\rho C}$ is the diffusion
coefficient, which has dimensions $[\alpha^2]= \frac{m^2}{s}$.

A similar line of reasoning for the heat flow in a
conduction plate leads to the two dimensional Heat Equation:

$$
\begin{equation}
\frac{\partial u}{\partial t}=\alpha^2\left(\frac{\partial^2 u}{\partial
x^2}+\frac{\partial^2 u}{\partial y^2}\right)
\end{equation}
$$(ref-intro-continued-6)

#### Fick's Law and Diffusion

Another way to arrive at the diffusion equation is to return to the
conservation law {eq}`ConservationLaw`, but instead of $u(x,t)$
representing the density of cars on a freeway, let us interpret $u$
as the concentration of molecules of a certain chemical in a stream
and $q$ the flux of these molecules. In this case the constitutive
relation between $q$ and $u$ is provided by Fick's law

$$
\begin{equation}
 q=-\alpha^2 \frac{\partial u}{\partial x}
\label{eq:Ficks_Law}
\end{equation}
$$(Ficks_Law)

Combining {eq}`Ficks_Law` with the conservation law
{eq}`ConservationLaw` we recover {eq}`Heat_equation`,
which in this context is known as the diffusion equation. That the
same equation aliases as the heat equation or the diffusion
equation stems from the distinct areas of application from which
these names arise.

#### The Drunkard's walk

Consider fruit flies having a density $u(x,t)$ at point $x$ at time
$t$ located on a row of trees that are spaced $\Delta x$ apart. We
assume that the fruit flies will jump to the tree to the left with a
probability $p$ and to the right with the same probability $p$. The
probability that the fruit files stay on the tree is $1-2p$.
Find an equation for the density of flies at $t+\Delta t$, i.e., $u(x,t+\Delta t)$.

```{figure} ../img/intro/drunkards_walk.png
:name: drunkards_walk
:alt: drunkards_walk
:align: center

Consider fruit flies having a density $u(x,t)$ located on a
row of trees that are spaced $\Delta x$ apart.
```

$$
\begin{eqnarray}
u(x,t+\Delta t)&=& pu(x+\Delta x,t)+(1-2p)u(x,t)+pu(x-\Delta x,t)\nonumber\\
&=& u(x,t)+p \Delta x\frac{[u(x+\Delta x,t)-u(x,t)]}{\Delta
x}-\frac{[u(x,t)-u(x-\Delta x,t)]}{\Delta x} \nonumber\\
&\simeq & u(x,t)+p\Delta x^2\frac{\left\{\frac{\partial u}{\partial
x}(x,t)-\frac{\partial u}{\partial x}(x-\Delta x,t)\right\}}{\Delta x}\nonumber\\
&\simeq & u(x,t)+p\Delta x^2\frac{\partial ^2u}{\partial x^2}\nonumber\\
\frac{u(x,t+\Delta t)-u(x,t)}{\Delta t}&\simeq & \left( p\frac{\Delta
x^2}{\Delta t}\right) \frac{\partial ^2u}{\partial x^2}\label{discrete_drunkard}
\end{eqnarray}
$$(discrete_drunkard)

Now choose the mesh and time sampling such that when taking the limit $\Delta x,\ \Delta t \rightarrow 0$ we obtain the limiting value $\left( p\frac{\Delta x^2}{\Delta t}\right)\rightarrow
\alpha^2$ so that {eq}`discrete_drunkard` reduces to

$$
\begin{equation}
 \frac{\partial u}{\partial t}=\alpha^2\frac{\partial^2 u}{\partial
x^2}\quad\mbox{The Heat Eq.} \label{eq:Diffusion_Equation}
\end{equation}
$$(ref-intro-continued-9)

Question: What is the Mean Absolute Deviation of a fruit fly?

```{figure} ../img/intro/DMeanDeviation.png
:name: DMeanDeviation.png
:alt: DMeanDeviation.png
:align: center

Consider the trajectory a single fruit fly in which it
takes steps of size $\pm\Delta x$ each time step $\Delta t$. What is the
mean deviation about its starting point?
```

$$
\begin{eqnarray*}
s_j&=&\pm\Delta x\\
t_j&=& j\Delta t
\end{eqnarray*}
$$(ref-intro-continued-10)

$$
\begin{eqnarray}
x_N&=&s_1+s_2+\cdots +s_N\sim 0\quad\mbox{Expected Value}\\
x_N^2&=& (s_1+\cdots + s_N)^2=s_1^2+\cdots +s_N^2+2(s_1s_2+\cdots
+s_{N-1}s_N)\nonumber\\
& &\hspace{1.1in}\sim N\Delta x^2\nonumber
\end{eqnarray}
$$(ref-intro-continued-11)

Therefore

$$
\begin{eqnarray*}
x_N^2&\simeq & \left(\frac{t_N}{\Delta t}\right) \Delta x^2=k^2t_N\\
|x_N|&\sim & k\sqrt{t_N}
\end{eqnarray*}
$$(ref-intro-continued-12)

```{figure} ../img/intro/drunkplot.png
:name: drunkplot
:alt: drunkplot
:align: center

Simulation with $N=1000$ trajectories for 200 steps and the
mean absolute deviation envelopes shown in red.
```
