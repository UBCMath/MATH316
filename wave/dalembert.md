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

# Interpretating D'Alembert's Solution in Space-Time

In this lecture we discuss the physical interpretation of the
D'Alembert solution in terms of space-time plots. In particular we
identify the left and right-moving characteristics as well as
the domain of dependence of a given point $(x_0,t_0)$ in space-time
and the region of influence of a given initial value specified at
the point ($x_1$,0). We discuss the evolution of a few simple pulses
and track the regions in space-time that are carved out by the
intersecting characteristics.

```{admonition} Key Concepts
The one dimensional Wave Equation;
D'Alembert's Solution, Characteristics, Domain of Dependence, Region
of Influence.
```

Reference Section: Boyce and Di Prima Section 10.7

## Space-Time Interpretation of D'Alembert's Solution

In this lecture we discuss the interpretation of D'Alembert's solution

$$
\begin{eqnarray}
u(x,t)=\frac{1}{2}\left[ u_0(x-ct)+u_0(x+ct)\right]
   +\frac{1}{2c}\int\limits_{x-ct}^{x+ct} v_0(s)ds  \label{DAlem}
\end{eqnarray}
$$(DAlem)

to the one dimensional wave equation

$$
\begin{eqnarray}
\frac{\pa^2 u}{\pa t^2}=c^2\frac{\pa^2 u}{\pa x^2} \label{eq1DW}
\end{eqnarray}
$$(ref1)

### Characteristics

In the $x-t$ plane the lines

$$
\begin{eqnarray}
x-ct=x_0\mbox{ and }x+ct=x_0
\end{eqnarray}
$$(ref2)

are called the characteristics that emanate from the point
$(x_0,0)$ in space-time (see figure {numref}{wave_characteristics}). Characteristics
are the lines (or curves of more general hyperbolic problems) along
which information is propagated by the equation. To interpret the
characteristic lines in the $x-t$ plane, it is useful to rewrite the
characteristic equations in the form

$$
\begin{eqnarray}\begin{array}{lcl}
x-ct=x_0 &\Rightarrow &t=\phantom{-}\dst\frac{1}{c}x-\frac{1}{c}x_0\\
\\
x+ct=x_0 &\Rightarrow &t=-\dst\frac{1}{c}x+\frac{1}{c}x_0\end{array}
\end{eqnarray}
$$(ref3)

```{figure} ../img/wave/wave_characteristics.png
:name: wave_characteristics
:alt: wave_characteristics
:align: center

The characteristics that emanate from.
```

### Region of Influence and Domain of Dependence

__Region of influence:__ The lines $x+ct=x_1$ and $x-ct=x_1$ bound the region of
influence of the function values at the initial point $(x_1,0)$.
Thus all the solution values $u(x,t)$ within this region can be
influenced by the value at the point  $(x_1,0)$.

__Domain of Dependence:__ The lines $x=x_0-ct_0$ and $x=x_0+ct_0$ that pass through
the point $(x_0,t_0)$ bound the domain of dependence. Thus the
solution $u(x_0,t_0)$ depends on all the function values in the
shaded region.

```{figure} ../img/wave/region_of_influence.png
:name: region_of_influence
:alt: region_of_influence
:align: center

Space-time Region of Influence of the point $(x_1,0)$ and
Domain of Dependence of the point $(x0,t0)$, both of which can be
determined from D'Alembert's Solution {eq}`DAlem`.
```

````{prf:example} A Rectangular Pulse

$$
\begin{eqnarray}
u(x,0)=\left\{\begin{array}{ll}1&|x|<1\\0&|x|>1\end{array}\right\}=u_0(x)
\label{RectanngularP}
\end{eqnarray}
$$(RectanngularP)

$$
\begin{eqnarray}
u(x,t)=\frac{1}{2}\left[ u_0(x-ct)+u_0(x+ct)\right]
\end{eqnarray}
$$(ref5)

Let $c=1$.

$\mathbf{t=\frac{1}{2}}$:

$$
\begin{eqnarray}\begin{array}{lclll}
x_r-\dst\frac{1}{2}=1 &\Rightarrow &x_r=\dst\frac{3}{2}&x_R+\dst\frac{1}{2}=1&x_R=\dst\frac{1}{2}\\
\\
x_\ell -\dst\frac{1}{2}=-1 &\Rightarrow &x_\ell =-\dst\frac{1}{2}
  &x_L+\dst\frac{1}{2}=-1 &x_L=-\dst\frac{3}{2}\end{array}
\end{eqnarray}
$$(ref6)

$\mathbf{t=1}$:

$$
\begin{eqnarray}\begin{array}{lcllcl}
x_r-1=1 &\Rightarrow &x_r=2&x_R+1=1&\Rightarrow &x_R=0\\
x_\ell -1=-1&\Rightarrow &x_\ell =0 &x_\ell +1=-1 &\Rightarrow
&x_L=-2\end{array}
\end{eqnarray}
$$(ref7)

$\mathbf{t=2}$:

$$
\begin{eqnarray}\begin{array}{lcllcl}
x_r-2=1 &\Rightarrow &x_r=3 &x_R+2=1 &\Rightarrow &x_R=-1\\
x_\ell -2=-1 &\Rightarrow &x_\ell =1 &x_L+2=-1 &\Rightarrow
&x_L=-3\end{array}
\end{eqnarray}
$$(ref8)
````

```{figure} ../img/wave/wave_spreader.png
:name: wave_spreader
:alt: wave_spreader
:align: center

Top: Space-time representation of the regions in which the
solution takes on different values for the rectangular pulse
{eq}`RectanngularP`. Bottom: Cross sections of the solution
$u(x,t)$ at times $t=0,\ 1/2c,\ 1/c, $, and $ t>1/c$
```
