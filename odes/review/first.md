# First Order Equations

## Separable Equations

A first order differential equation is **separable** if it is of the form 

$$
\frac{dy}{dx} = P(x)Q(y)
$$(ref-odes-review-first-0)

Solve a separable equation by dividing by $Q(y)$ and integrating with respect to $x$

$$
\int\frac{dy}{Q(y)} = \int P(x) dx
$$(ref-odes-review-first-1)

This is called the method of **separation of variables**. For example, find the general solution of the equation

$$
\frac{dy}{dx} = xy
$$(ref-odes-review-first-2)

Apply the method of separation of variables

$$
\int \frac{dy}{y} &= \int x \, dx \\
\ln|y| &= \frac{x}{2} + C \\
y &= A e^{x/2}
$$(ref-odes-review-first-3)

where $A = e^C \in \mathbb{R}$. We found an explicit closed form solution $y(x)$ in this case however this is not always possible. For example, find the general solution of the equation

$$
\frac{dy}{dx} = \frac{4y}{x(y-3)}
$$(ref-odes-review-first-4)

by applying the method of separation of variables

$$
\int \left( 1 - \frac{3}{y} \right) \, dy &= 4 \int \frac{dx}{x} \\
y - 3 \ln|y| &= 4 \ln|x| + C \\
y &= \ln (x^4 y^3) + C \\
A x^4 y^3 &= e^y
$$(ref-odes-review-first-5)

where $A = e^C \in \mathbb{R}$. We cannot solve for $y$ therefore this defines an implicit solution.

## Linear Equations

A first order differential equation is **linear** if it is of the form

$$
y'(x) + P(x)y = Q(x)
$$(ref-odes-review-first-6)

Can we find a function $F(x)$ to multiply {eq}`ref-odes-review-first-6` by in order to turn the left hand side into a derivative of a product

$$
F(x)y' + F(x)P(x)y = F(x)Q(x)
$$(ref-odes-review-first-7)

$$
(F(x)y)' = F(x)y' + F'(x)y = F(x)Q(x)
$$(ref-odes-review-first-8)

So let $F'(x) = F(x)P(x)$ which is a separable equation and so we apply separable of variables

$$
\int \frac{dF}{F} &= \int P(x)dx \\
\ln F &= \int P(x)dx \\
F(x) &= e^{\int P(x)dx}
$$(ref-odes-review-first-9)

The function $F(x)$ is called **integrating factor**

$$
F(x) = e^{\int P(x) dx}
$$(ref-odes-review-first-10)

This leads us to an expression for the general solution of {eq}`ref-odes-review-first-6`

$$
F(x) y' + F(x)P(x)y &= F(x)Q(x) \\
\left( F(x) y \right)' &= F(x) Q(x) \\
y(x) &= \frac{1}{F(x)} \left( \int F(x)Q(x)dx + C \right)
$$(ref-odes-review-first-11)

For example, let's apply the integrating factor to find the general solution of the linear equation

$$
y' + 2y = 0
$$(ref-odes-review-first-12)

The integrating factor is

$$
F(x) = e^{2x}
$$(ref-odes-review-first-13)

Therefore we find the solution

$$
e^{2x}y' + 2e^{2x}y &= 0 \\
\left( e^{2x}y \right)' &= 0 \\
e^{2x}y &= C \\
y(x) &= C e^{-2x}
$$(ref-odes-review-first-14)

Let's try another example. Find the unique solution of the equation

$$
\frac{dy}{dx} + \cot(x) y = 5 e^{\cos(x)} \ \ , \ \ \ y(\pi/2) = -4
$$(ref-odes-review-first-15)

Identify the functions $P(x)$ and $Q(x)$ and find the integrating factor

$$
P(x) &= \cot(x) \\
Q(x) &= 5 e^{\cos(x)} \\
F(x) &= e^{\int \cot(x)dx} = e^{\ln(\sin(x))} = \sin(x)
$$(ref-odes-review-first-16)

Compute the general solution

$$
\sin(x) y' + \cos(x) y &= 5 e^{\cos(x)}\sin(x) \\
(\sin(x) y)' &= 5 e^{\cos(x)}\sin(x) \\
\sin(x) y &= -5 e^{\cos(x)} + C\\
y(x) &= - \frac{5e^{\cos(x)} - C}{\sin(x)} \\
$$

Plug in the initial value to find $C$

$$
-4 &= y(\pi /2)= - \frac{5 - C}{1} \ \ \Rightarrow \ \ C = 1 \\
$$(ref-odes-review-first-17)

The unique solution of {eq}`ref-odes-review-first-15` is

$$
y(x) = \frac{1-5 e^{\cos(x)}}{\sin(x)}
$$(ref-odes-review-first-18)