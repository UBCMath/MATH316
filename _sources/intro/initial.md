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

# Introduction To Partial Differential Equations

In this lecture we will introduce the three basic partial
differential equations we consider in this course. After briefly
discussing the classification of these equations we go through the
modeling process to arrive at these three equations from a number of
different physical situations. We deliberately explore the different
paths to arrive at the same partial differential equations to
emphasize way in which the models from disparate applications share
the same features.

```{admonition} Key Concepts
Partial Differential Equations (PDEs);
Elliptic, Parabolic, Hyperbolic PDEs; The heat Equation, The Wave
Equation, and Laplace's Equation, Modeling and Derivation of PDEs.
```

## Introduction To PDEs

### Classification Of PDEs

___Ordinary Differential Equations (ODE):___

Equations which define
functions of a single independent variable by prescribing a
relationship between the values of the function and its derivatives.

````{prf:example}
:label: example-intro-initial-0 A Nonlinear Second-Order ODE
$$
\begin{equation}
y''(x)+e^{y(x)}=0.
\end{equation}
$$(ref-intro-initial-0)
````

___Partial Differential Equations (PDEs):___

Involve multivariable functions $u(x,t)$, $u(x,y)$ that are determined by
prescribing a relationship between the function value and its
partial derivatives.

````{prf:example}
:label: example-intro-initial-1 A Linear First-Order PDE
$$
\begin{equation}a(x,y) \frac{\partial }{\partial
x}u(x,y)+b(x,y) \frac{\partial }{\partial y}u(x,y)=c(x,y)
\end{equation}
$$(ref-intro-initial-1)
````

````{prf:example}
:label: example-intro-initial-2 A Nonlinear First-Order PDE
$$
\begin{equation}a(x,y,u) \frac{\partial }{\partial
x}u(x,y)+b(x,y,u) \frac{\partial }{\partial y}u(x,y)=c(x,y,u)
\end{equation}
$$(ref-intro-initial-2)
````

#### Some Classic Linear Second-Order PDEs

$$
{\small
\begin{array}{|c|c|c|c|}
\hline &&&\\
\mbox{Quadric} &\mbox{Classification} &\mbox{Eq.} &\mbox{Name}\\
&&&\\
\hline &&&\\
T=X^2&\mbox{Parabolic}&\frac{\partial u(x,t)}{\partial
t}=\frac{\partial^2 u(x,t) }{\partial x^2}&\mbox{The Heat Equation or Diffusion Equation}\\
&&&\\
\hline &&&\\
X^2+Y^2=k&\mbox{Elliptic}&\frac{\partial^2 u(x,y) }{\partial
x^2}+\frac{\partial^2 u(x,y)}{\partial y^2}=f(x,y)
&\begin{array}{ll}
\mbox{Poisson's Equation}\, f\not\equiv 0\\
\mbox{Laplace's Equation}\, f=0\end{array}\\
&&&\\
\hline &&&\\
T^2-c^2X^2=k&\mbox{Hyperbolic}&\frac{\partial^2 u(x,t)}{\partial
t^2}-c^2\frac{\partial^2 u(x,t)}{\partial x^2}=0&\mbox{The Wave Equation}\\
&&&\\
\hline
\end{array} }
$$

By analogy with quadric surfaces $aX^2+2bXY+c^2Y^2+\cdots =k$
that can be reduced to a standard form by coordinate rotation, the
most general linear 2nd order PDE

$$
\begin{equation}au_{xx}+2bu_{xy}+cu_{yy}+\cdots\end{equation}
$$(ref-intro-initial-3)

can be reduced by a transformation of coordinates to one of the
Heat, Laplace or the Wave Equations.

### Modeling And The Derivation Of PDE

Mathematical modeling is the process of writing down an equation or
a system of equations that describe a particular physical, chemical,
biological, or economic system that we wish to understand at a more
fundamental level and whose behavior we would like to predict and
perhaps even control. Mathematical modeling is an art-form in which
a loose set of tools are applied to arrive at a self-consistent
model, which can give a faithful representation of the behavior of
the target system - sometimes with startling results. Tools that are
typically used to build mathematical models include: conservation
principles and balance laws that must obviously govern the behavior
of the system we are trying to describe, e.g. conservation of mass,
fluid, chemicals, fruit flies or balance of linear momentum,
Newton's Second and Third Laws of Mechanics. 

In this process it is
very important to be very mindful of the dimensions of the variables
that we define so that we do not commit the cardinal sin of "adding
apples to oranges". Dimensional analysis, rather than being a mere
check for consistency, has evolved to an extremely useful sub-field
of ODE and PDE that can be used to derive properties of certain
solutions and even to identify special solutions that could not be
obtained by other techniques. Other reality checks in the modeling
process are obtained by ensuring that the model does not violate
some very fundamental physical or economic principle - such as the
second law of thermodynamics or the postulate of a liquid market.

Modeling is an extremely broad topic, which arguably includes all of physics,
physical chemistry, mathematical biology, and about which many books
have been written. Therefore we will not have time to explore this
topic in much detail in this course. We will, however, explore a few
examples to illustrate the modeling process. One aspect of the
remainder of this lecture to which you should pay particular
attention is the way in which we can arrive at precisely the same
equation in spite of the fact that we are considering completely
different physical systems with very different meanings attached to
the dependent and independent variables. Thus modeling is a
tremendously unifying process, which can highlight the fundamental
similarities between the behavior of seemingly disparate physical
systems. In fact, we will only be studying three equations in this
course! However, the richness of the diverse applications of these
few equations is what makes Applied Mathematics so interesting.
Indeed, it is the reason that Applied Mathematicians are in such
high demand in almost every field of industry from: geoengineering,
e.g., mining, extracting petroleum, geophysical prospecting; to
every branch of engineering, e.g. to design more efficient circuits,
medical devices and imaging techniques, or to designing safer
aircraft. The focus of this course will be on what comes after the
model has been built, i.e., given a mathematical model how do we
find a solution? Given this emphasis you may be tempted to forget
about the modeling aspect of the PDE you will find that you can
derive much insight about the behavior of the solutions by keeping
in mind the physical meaning behind the variables. For example, a
simple ``thought experiment", with a mental picture of one of the
physical systems to which a given PDE applies, can be used to check
a solution that you have derived to see if it makes sense.

### A One Dimensional Conservation Law

#### Traffic Flow on a Highway

Consider the traffic flow on a highway and let $u(x,t)$ be the 
_density_ of cars at $x$ at time $t$.

$$
\begin{equation}
[u]=\, \mbox{$\# $ of cars/unit length.}
\end{equation}
$$(ref-intro-initial-4)

Let $q(x,t)$ be the flux of cars at $x$ at time $t$.

$$
\begin{equation}
[q]=\, \mbox{$\# $ of cars/unit time.}
\end{equation}
$$(ref-intro-initial-5)

```{figure} ../img/intro/traffic_flow_conservation_law.png
:name: traffic_flow_conservation_law
:alt: traffic_flow_conservation_law
:align: center

Traffic flow along the $x$ axis with density $u(x,t)$
(cars/unit length) and flux $q(x,t)$ (cars/second) at $x$ \& instant
$t$.
```

Now the change in the number of cars within the interval
$[x,x+\Delta x]$ is approximated by

$$
\begin{equation}
\Delta u \Delta x = \{ u(x,t+\Delta t)-u(x,t)\}\Delta x\simeq\{
q(x,t)-q(x+\Delta x,t)\}\Delta t  \label{eq:changeinCarsEQ}
\end{equation}
$$(changeinCarsEQ)

Now divide {eq}`changeinCarsEQ` by $\Delta t \Delta x$ and let
$\Delta x\rightarrow 0$ and $\Delta t\rightarrow 0$ and we obtain the following
conservation law PDE:

$$
\begin{equation}
 \frac{\partial u}{\partial t}+\frac{\partial q}{\partial x}=0
\label{eq:ConservationLaw}
\end{equation}
$$(ConservationLaw)

 This limiting process is frequently referred to as ''taking the continuum
limit".

````{prf:observation}
:label: observation-intro-initial-0
This partial differential equation represents the
conservation of a quantity $u(x,y)$ that is subject to a flux
$q(x,t)$, which is why it is called a conservation law. Depending on
the context and the definitions of $u$ and $q$, the conservation law
can be used to represent the following physical situations (among
many) in which quantities are conserved:

- conservation of cars

- conservation of heat

- conservation of chemicals.

- conservation of fluid.

We observe that the conservation law relates the gradients
of two distinct quantities $u$ and $q$. In order to have enough
information to solve for one of the variables, we need to provide
another equation. This is sometimes provided by what is known as
an equation of state in thermodynamics or a constitutive
relation in continuum mechanics. For example, how does $q$ change
with $u$ or its derivatives, i.e., $q=q(u)$, $q=q(x,t,u)$, or
$q=q(x,t,u,\frac{\partial u}{\partial x})$?

```{figure} ../img/intro/galilean_transformation.png
:name: galilean_transformation
:alt: galilean_transformation
:align: center

The Galilean transformation of coordinates from $x$ to
$x'=x-ct $
```
````

#### Application: Convection and the First Order Wave Equation

Assume that the flux of cars $q$ increases linearly with the density
of cars $u$, i.e., $q=cu,\ c>0$, then it follows that

$$
\begin{equation}
 \frac{\partial u}{\partial t}+c\frac{\partial u}{\partial x}=0
\label{eq:RightMovingWave}
\end{equation}
$$(RightMovingWave)

But this is just a wave equation. To see this consider
two coordinate systems $Ox$ and $O'x'$. Assume that at time $t=0$
the coordinate systems $Ox$ and $O'x'$ are coincident and that
$O'x'$ moves at a speed $c$ relative to coordinate system $Ox$ and
directed toward increasing $x$. To make the situation more realistic
assume that the moving coordinate system is attached to a wave whose
shape is shown in figure
{numref}`galilean_transformation` and that the
red "surfer" is riding with the wave. We assume that there is a
blue "observer" attached to the fixed coordinate system $Ox$. At
time $t=0$, since the two coordinate systems $Ox$ and $O'x'$ were
coincident, the blue and red observers were at the same place.
According to the surfer the functional form of the wave represented
by the function $f(x')$ stays the same throughout the motion. We
observe that the distance between the $O'$ and the vertical red line
remains $x'$ throughout the motion, while the distance $x$ from the
centre $O$ of the stationary coordinate system is related to $x'$ by
the so-called Galilean Transformation:

$$
\begin{equation}
 x'=x-ct
\label{eq:Galilean_Transformation}
\end{equation}
$$(Galilean_Transformation)

Thus according to the stationary observer, the functional
form of the wave varies in space-time according to:

$$
\begin{equation}
 f(x')=f(x-ct)
\label{eq:Stationary_Observer}
\end{equation}
$$(Stationary_Observer)

Motivated by this property of a wave/signal moving to the
right at a constant speed $c$, we are led to consider the following
guess for a solution to {eq}`RightMovingWave`:

$$
\begin{eqnarray}
\mbox{We guess that}\quad u(x,t)&=&f(x-ct)\quad\mbox{solves}\quad u_t+cu_x=0 \label{sol:RightMovingWave}\\
\mbox{Take derivatives }\quad u_t&=&-cf'\quad u_x=f' \nonumber\\
\mbox{Therefore}\quad u_t+cu_x&=&-cf'+cf'=0, \quad \mbox{which
implies}\quad u(x,t)=f(x-ct)\quad \mbox{solves
{eq}`RightMovingWave`}.\nonumber
\end{eqnarray}
$$(RightMovingWaveSol)

Thus $u_t+cu_x=0$ has solutions of the form $u(x,t)=f(x-ct)$ for any
sufficiently differentiable $f$, each of which represents a right
moving wave of a given shape.

````{prf:observation}
:label: observation-intro-initial-1

- _A judicious guess:_ Because {eq}`RightMovingWave`
comprises a linear combination of a time derivative $\frac{\partial }{\partial
t}$ and a spatial derivative $\frac{\partial }{\partial x}$, we might expect
to find a solution of the form of an exponential of a linear
function of $x$ and $t$, since either derivative of such a function
is in the form of a constant times the exponential. We therefore
consider the trial solution of the form

$$
\begin{equation}
u(x,t) = \{\rm\ e\}^{\{\rm\ i\} kx + \sigma t}
 \label{eq:ExponentialSpaceTimeSolution}
\end{equation}
$$(ExponentialSpaceTimeSolution)

Substituting {eq}`ExponentialSpaceTimeSolution` into
{eq}`RightMovingWave` we obtain

$$
\begin{equation}
\left(\frac{\partial}{\partial t}+c\frac{\partial}{\partial x}\right) \{\rm\ e\}^{\{\rm\ i\} kx +
\sigma t} = \{\sigma + \{\rm\ i\} k c\}\{\rm\ e\}^{\{\rm\ i\} kx + \sigma t},
 \label{eq:1DWaveEqDispersionRelationDerivation}
\end{equation}
$$(ref-intro-initial-13)

which is a solution of provided $\sigma$ and $k$ satisfy the
following "dispersion relation"

$$
\begin{equation}
\sigma = - \{\rm\ i\} k c
 \label{eq:1DWaveEqDispersionRelation}
\end{equation}
$$(ref-intro-initial-14)

Thus we have obtained a special case of the solution derived in
{eq}`RightMovingWaveSol`

$$
\begin{equation}
u(x,t) = g(x-ct) = \{\rm\ e\}^{\{\rm\ i\} k \left(x - c t \right)}
 \label{eq:ExponentialSpaceTimeSolution_g}
\end{equation}
$$(ExponentialSpaceTimeSolution_g)

- _Shocking - a nonlinear wave equation:_ What happens if
$q$, instead of increasing linearly with $u$, can behave in a
nonlinear way, i.e.,

$$
\begin{equation}
q(u) = h(u), \mbox{ where $h$ is some given function of u}
 \label{eq:NonlinearWaveEqFlux}
\end{equation}
$$(ref-intro-initial-16)

Combining this with the conservation law
{eq}`ConservationLaw` can be written in the form

$$
\begin{equation}
 \frac{\partial u}{\partial t}+h'(u) \frac{\partial u}{\partial x}=0
\label{eq:NonlinearWaveEq}
\end{equation}
$$(ref-intro-initial-17)
We note that since the wave speed $c~h'(u)$ can vary in space, it is
possible for certain initial conditions and functions $h$ to have the
waves that initiate for more negative values of $x$ to crash into
waves that initiate for more positive values of $x$. This phenomenon
is resolved in wave mechanics by the formation of a shock
wave, which represents a special solution to this over-specified
situation in which there are potentially multiple values of the
solution at certain points in the domain. This is the same
phenomenon that occurs with formation of supersonic shock waves by
aircraft or by the cracking of a whip.

- _A left moving wave:_ What happens if $q=-cu$? In this case

$$
\begin{equation}
 \frac{\partial u}{\partial t}-c\frac{\partial u}{\partial x}=0
\label{eq:LeftMovingWave}
\end{equation}
$$(LeftMovingWave)

We leave it as an exercise to show, in a similar way to the
procedure used for the right moving wave, that
{eq}`LeftMovingWave` has a solution that represents a left
moving wave.

- _The second order wave equation:_ Note that if we apply the left 
$\frac{\partial }{\partial t}-c\frac{\partial }{\partial x}$ and right 
$\frac{\partial }{\partial t}+c\frac{\partial }{\partial x}$ moving wave
operators in succession, we obtain

$$
\begin{equation}
\left(\frac{\partial}{\partial t}+c\frac{\partial}{\partial
x}\right)\left(\frac{\partial}{\partial t}-c\frac{\partial}{\partial x}\right)
u(x,t)=\frac{\partial^2u}{\partial t^2} -c^2\frac{\partial^2 u}{\partial x^2}=0,
\label{2nd_oder_wave_eq}
\end{equation}
$$(2nd_oder_wave_eq)

which is the second order wave equation that has both left and right
moving wave solutions (we will return to this later in the context
of acoustic waves in a solid bar).
````

#### Application: the Convection-Diffusion equation

Consider the traffic flowing down the highway as shown in figure
{numref}`traffic_flow_conservation_law` and assume that the flux $q$
increases linearly with the car density $u$. Now assume some agency
from the driver in which she responds to an increase in the density
of traffic by decreasing her speed, which results in a decrease in
the flux of cars locally. This situation can be represented by a
flux function of the form

$$
\begin{equation}
 q=cu-D \frac{\partial u}{\partial x}
\label{eq:ConvectionDiffusionFlux}
\end{equation}
$$(ConvectionDiffusionFlux)

Combining {eq}`ConvectionDiffusionFlux` with
{eq}`ConservationLaw` we obtain the convection-diffusion
equation

$$
\begin{equation}
 \frac{\partial u}{\partial t}+c \frac{\partial u}{\partial x}=D \frac{\partial^2 u}{\partial x^2}
\label{eq:Convection_Diffusion_Eq}
\end{equation}
$$(Convection_Diffusion_Eq)

````{prf:observation}
:label: observation-intro-initial-2
- A second order parabolic PDE: Considering the highest
derivatives that appear in each of the independent variables $x$ and
$t$ we observe (from the chart at the beginning of this lecture)
that the convection-diffusion is classified as a parabolic PDE.

- A moving coordinate system: Introduce the transformation $u(x,t)=U(\xi ,t),$ where $\xi =x-ct$ and
show that $U(\xi ,t)$ satisfies the diffusion equation
$U_{t}=DU_{\xi \xi }.$ Can you interpret the removal of the
convection term physically?

- The dispersion relation and stability: Consider solutions of {eq}`Convection_Diffusion_Eq` of the form $u(x,t)=\mathrm{
e}^{\mathrm{i}kx+\sigma t}.$ Determine the associated dispersion
relation $\sigma =\sigma (k)$. Using the dispersion relation
determine if the solution stable when $D>0$ or when $D<0?$
````
