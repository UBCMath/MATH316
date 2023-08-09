# Cauchy-Euler Equations

## Substitution $x = e^t$

A second order linear equation of the form

$$
L[y] = x^2y'' + \alpha xy' + \beta y = 0 \ \ , \ \ \ \alpha , \beta \in \mathbb{R}
$$(ref-odes-review-cauchy-0)

is called a **Cauchy-Euler equation** (or an **equidimensional equation**). We can solve a Cauchy-Euler equation by using a **substitution** to transform it into a second order equation with constant coefficients. In particular, let

$$
x &= e^t \\
t &= \ln x
$$(ref-odes-review-cauchy-1)

The chain rule gives the relation

$$
\frac{d}{dx} = \frac{d}{dt} \frac{dt}{dx} \ \ \Rightarrow \ \ \frac{d}{dt} = x\frac{d}{dx}
$$(ref-odes-review-cauchy-2)

and also

$$
\frac{d^2}{dt^2} &= x \frac{d}{dx} \left( x \frac{d}{dx} \right) \\
&= x^2 \frac{d^2}{dx^2} + x\frac{d}{dx}
$$(ref-odes-review-cauchy-3)

In particular, we have

$$
x^2 \frac{d^2}{dx^2} = \frac{d^2}{dt^2} - \frac{d}{dt}
$$(ref-odes-review-cauchy-4)

Therefore, using the substitution {eq}`ref-odes-review-cauchy-1`, the equation {eq}`ref-odes-review-cauchy-0` becomes a second order equation with constant coefficients

$$
\frac{d^2y}{dt^2} - \frac{dy}{dt} + \alpha \frac{dy}{dt} + \beta y &= 0 \\
\frac{d^2y}{dt^2} + (\alpha - 1) \frac{dy}{dt} + \beta y &= 0
$$(ref-odes-review-cauchy-5)

We know that solutions of {eq}`ref-odes-review-cauchy-5` are of the form $y = e^{rt}$ therefore to solve {eq}`ref-odes-review-cauchy-0` we use the substitution {eq}`ref-odes-review-cauchy-1` and let

$$
y &= x^r \\
y' &= rx^{r-1} \\
y'' &= r(r-1)x^{r-2}
$$(ref-odes-review-cauchy-6)

and {eq}`ref-odes-review-cauchy-0` becomes

$$
\left( r(r-1) + \alpha r + \beta\right) x^r = 0
$$(ref-odes-review-cauchy-7)

The **characteristic polynomial** of the Cauchy-Euler equation {eq}`ref-odes-review-cauchy-0` is

$$
r^2 + (\alpha - 1)r + \beta = 0
$$(ref-odes-review-cauchy-8)

with roots

$$
r_1 &= \frac{1 - \alpha + \sqrt{(\alpha - 1)^2 - 4\beta}}{2} \\
r_2 &= \frac{1 - \alpha - \sqrt{(\alpha - 1)^2 - 4\beta}}{2}
$$(ref-odes-review-cauchy-9)

We know that the form of the solution of {eq}`ref-odes-review-cauchy-0` depends on the **discriminant**

$$
\Delta = (\alpha - 1)^2 - 4\beta
$$(ref-odes-review-cauchy-10)

## Case I: Real Distinct Roots

If $\Delta > 0$ then the roots $r_1$ and $r_2$  are real and distinct and the general solution of {eq}`ref-odes-review-cauchy-0` is

$$
y(x) = c_1x^{r_1} + c_2x^{r_2} \ , \ \ c_1,c_2 \in \mathbb{R}
$$(ref-odes-review-cauchy-11)

Note that if either $r_1 < 0$ or $r_2 < 0$ then $|y| \to \infty$ as $x \to 0$.

## Case II: Real Repeated Roots

If $\Delta = 0$ then the roots are real and repeated and 

$$
r_1 = r_2 = \frac{1 - \alpha}{2}
$$(ref-odes-review-cauchy-12)

We obtain only one solution in this case

$$
y = c_1 x^{r_1} \ , \ \ c_1 \in \mathbb{R}
$$(ref-odes-review-cauchy-13)

To get a second solution we use the method introduced in {ref}`ref-odes-review-second` (see {ref}`ref-odes-review-second-repeated`) in which we differentiate with respect to the parameter $r$

$$
\frac{\partial}{\partial r} \left( x^r \right) = x^r \ln x
$$(ref-odes-review-cauchy-14)

Let's verify that $y(x) = x^{r_1} \ln x$ is a solution of {eq}`ref-odes-review-cauchy-0`. Compute the derivatives with respect to $x$

$$
(x^r \ln x)' &= rx^{r-1}\ln x + x^{r-1} \\
(x^r \ln x)'' &= r(r-1)x^r\ln x + rx^{r-2} + (r-1)x^{r-2}
$$(ref-odes-review-cauchy-15)

and plug in to {eq}`ref-odes-review-cauchy-0` to find

$$
L[x^r \ln x] &= x^2 (x^r \ln x)'' + \alpha x(x^r \ln x)' + \beta (x^r \ln x) \\
&= x^2 \left( r(r-1)x^r\ln x + rx^{r-2} + (r-1)x^{r-2} \right) \\
& \quad \quad + \ \alpha x \left( rx^{r-1}\ln x + x^{r-1} \right) + \beta (x^r \ln x) \\
&= \left( r^2 + (\alpha - 1)r + \beta \right) x^r \ln x + \left( 2r - 1 + \alpha \right) x^r
$$(ref-odes-review-cauchy-16)

Therefore $x^{r_1} \ln x$ is indeed a solution. Finally, the general solution of {eq}`ref-odes-review-cauchy-0` when $\Delta = 0$ is

$$
y(x)=(c_1 + c_2\ln x)x^{r_1} \ , \ \ c_1,c_2 \in \mathbb{R}
$$(ref-odes-review-cauchy-17)

## Case III: Complex Conjugate Roots

If $\Delta < 0$ then then $r_1$ and $r_2$ complex conjugate roots. In particular, write $r_1 = \lambda + i \mu$ and $r_2 = \lambda - i \mu$ where

$$
\lambda = \frac{1 - \alpha}{2} \hspace{20mm} \mu = \frac{ \sqrt{4 \beta - (\alpha - 1)^2}}{2}
$$(ref-odes-review-cauchy-18)

The general solution of {eq}`ref-odes-review-cauchy-0` when $\Delta < 0$ is given by

$$
y(x) &= c_1 x^{(\lambda + i\mu)} + c_2x^{(\lambda - i\mu)} \\
&= c_1 e^{(\lambda + i\mu) \ln x} + c_2 e^{(\lambda - i\mu) \ln x} \\
&= x^{\lambda} \left( c_1 e^{i \mu \ln x} + c_2 e^{-i\mu\ln x} \right)
$$(ref-odes-review-cauchy-19)

and fincally

$$
y(x) = A_1 x^\lambda \cos(\mu \ln x) + A_2 x^{\lambda} \sin(\mu \ln x) \ \ , \ \ \ A_1,A_2 \in \mathbb{R}
$$(ref-odes-review-cauchy-20)

If we consider the range $x < 0$ then we replace $x$ in the general solution {eq}`ref-odes-review-cauchy-20` by $|x|$.

Apply the Wronskian test to show that the two solutions in {eq}`ref-odes-review-cauchy-20` are linearly independent. The Wronskian is

$$
w(y_1,y_2) = \left| \begin{array}{cc} y_1 &y_2 \\ y_1' & y_2' \end{array} \right| = y_1 y_2' - y_1' y_2
$$(ref-odes-review-cauchy-21)

and the functions $y_1$ and $y_2$ are linearly independent if $w(y_1,y_2) \not = 0$. Let $y_1 = x^{\lambda} \cos (\mu\ln x)$ and $y_2 = x^{\lambda} \sin (\mu\ln x)$ and compute

$$
w(y_1,y_2) &= x^{\lambda} \cos(\mu \ln x) \left( x^{\lambda} \ln x \sin(\mu\ln x) + x^{\lambda - 1} \cos(\mu \ln x) \mu \right) \\
& \quad \quad - \ x^{\lambda} \sin(\mu \ln x) \left( x^{\lambda} \ln x \cos(\mu \ln x) - x^{\lambda - 1} \sin(\mu \ln x) \mu \right) \\
&= \mu x^{2\lambda - 1}
$$(ref-odes-review-cauchy-22)

Therefore we conclude the solutions are linearly independent.

## Examples

**Example 1.** Find the unique solution of the equation

$$
x^2 y'' - xy' - 2y = 0 \ , \ \ y(1) = 0 \ , \ \ y'(1) = 1
$$(ref-odes-review-cauchy-23)

The characteristic polynomial is $r^2 - 2r - 2 = 0$ and the roots are

$$
r_1 = 1 + \sqrt{3} \hspace{20mm} r_2 = 1 - \sqrt{3} 
$$(ref-odes-review-cauchy-24)

therefore the generral solution of {eq}`ref-odes-review-cauchy-23` is

$$
y(x) = c_1 x^{1 + \sqrt{3}} + c_2 x^{1 - \sqrt{3}}
$$(ref-odes-review-cauchy-25)

The initial condition $y(1)=0$ implies $y(1) = c_1 + c_2 = 0$ therfore $c_2 = -c_1$

$$
y(x) = c_1 \left( x^{1 + \sqrt{3}} - x^{1 - \sqrt{3}} \right)
$$(ref-odes-review-cauchy-26)

Compute the derivative

$$
y'(x) = c_1 \left( \left( 1 + \sqrt{3} \right) x^{\sqrt{3}} - \left( 1 - \sqrt{3} \right) x^{-\sqrt{3}} \right)
$$(ref-odes-review-cauchy-27)

and the condition $y'(1) = 1$ yields $c_1 = 1/(2 \sqrt{3})$. The unique solution of {eq}`ref-odes-review-cauchy-23` is given by

$$
y(x) = \frac{1}{2\sqrt{3}} \left(x^{1 + \sqrt{3}} - x^{1 - \sqrt{3}} \right)
$$(ref-odes-review-cauchy-28)

**Example 2.** Find the unique solution of the equation

$$
x^2 y'' - 3xy' + 4y = 0 \ , \ \ y(1) = 1 \ , \ \ y'(1) = 0
$$(ref-odes-review-cauchy-29)

The characteristic polynomial is $r^2 - 4r + 4 = 0$ and the roots are repeated $r_1 = r_2 = 2$ therefore the general solution is

$$
y(x) = c_1 x^2 + c_2 x^2 \ln x
$$(ref-odes-review-cauchy-30)

The initial condition $y(1) = 1$ implies $c_1 = 1$. Compute the derivative

$$
y'(x) = 2x + c_2 \left( 2x \ln x + x \right)
$$(ref-odes-review-cauchy-31)

and the condition $y'(1) = 0$ yields $c_2 = -2$. The unique solution of {eq}`ref-odes-review-cauchy-29` is given by

$$
y(x) = x^2 \left( 1 - 2 \ln x \right)
$$(ref-odes-review-cauchy-32)