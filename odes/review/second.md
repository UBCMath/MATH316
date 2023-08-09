(ref-odes-review-second)=
# Second Order Constant Coefficient Equations

## Characteristic Polynomial

A **second order constant coefficient homogeneous equation** is of the form

$$
ay'' + by' + cy = 0 \ \ , \ \ \ a,b,c \in \mathbb{R}
$$(ref-odes-review-second-0)

Introduce notation $L$ for the linear operator

$$
L[y] = ay'' + by' + cy 
$$(ref-odes-review-second-1)

Let $y(x) = e^{rx}$ and compute

$$
L \left[ e^{rx} \right] = (ar^2 + br + c) e^{rx}
$$(ref-odes-review-second-2)

Therefore $L[y] = 0$ provided $ar^2 + br + c = 0$ and solutions of {eq}`ref-odes-review-second-0` are determined by roots of the **characteristic polynomial**

$$
g(r) = ar^2 + br + c = 0
$$(ref-odes-review-second-3)

Denote the roots by

$$
r_1 = \frac{- b + \sqrt{b^2 - 4ac}}{2a} \ \ \ \text{and} \ \ \ r_2 = \frac{- b - \sqrt{b^2 - 4ac}}{2a}
$$(ref-odes-review-second-4)

and let

$$
\Delta = b^2 - 4ac
$$(ref-odes-review-second-5)

be the **discriminant** of the characteristic polynomial.

## Case I: Real Distinct Roots

If $\Delta > 0$ then $r_1$ and $r_2$ are **real distinct roots**: $r_1,r_2 \in \mathbb{R}$ and $r_1 \ne r_2$. The general solution of {eq}`ref-odes-review-second-0` in this case is

$$
y(x)=c_1 e^{r_1 x} + c_2 e^{r_2x} \ \ , \ \ \ c_1 , c_2 \in \mathbb{R}
$$(ref-odes-review-second-6)

(ref-odes-review-second-repeated)=
## Case II: Real Repeated Roots

If $\Delta = 0$ then $r_1$ and $r_2$ are **real repeated roots**: $r_1,r_2 \in \mathbb{R}$ and $r_1 = r_2$. In this case, we obtain only one solution

$$
y(x) = e^{r_1 x}
$$(ref-odes-review-second-7)

How do we get a second solution? We'll consider two methods to prove the existence of an another solution.

**First Method: Perturbation of the repeated root**

Consider a small perturbation $\epsilon$ to the repeated root (see {numref}`perturbed_roots`) and the perturbed characteristic polynomial

$$
g(r) = a(r - (r_1 - \epsilon))(r - (r_1 + \epsilon)) = a \left( (r - r_1)^2 - \epsilon^2 \right)
$$(ref-odes-review-second-8)

In this case, the two very close but distinct roots of $g(r) = 0$ are given by $r_1 + \epsilon$ and $r_1 - \epsilon$ and the general solution of {eq}`ref-odes-review-second-0` is

$$
y(x)= c_1 e^{(r_1 + \epsilon)x} + c_2 e^{(r_1 - \epsilon)x} \ \ , \ \ \ c_1 , c_2 \in \mathbb{R}
$$(ref-odes-review-second-9)

Choose a special solution by selecting

$$
c_1 = \frac{1}{2\epsilon} = -c_2
$$(ref-odes-review-second-10)

and we obtain a family of solutions that depend on the small parameter $\epsilon$ (see {numref}`perturbed_roots`)

$$
y(x,\epsilon) = \frac{e^{(r_1 + \epsilon)x} - e^{(r_1 - \epsilon)x}}{2\epsilon}
$$(ref-odes-review-second-11)

Take the limit as $\epsilon \rightarrow 0$ using L'Hospital's Rule to find the limiting solution

$$
y(x,\epsilon) = e^{r_1x} \left( \frac{e^{\epsilon x} - e^{-\epsilon x}}{2\epsilon} \right) \longrightarrow x e^{r_1x}
\hspace{10mm} \text{as} \ \ \epsilon \rightarrow 0
$$(ref-odes-review-second-12)

Therefore we have found that

$$
y(x) = x e^{r_1 x}
$$(ref-odes-review-second-13)

is another solution of {eq}`ref-odes-review-second-0` in the case $\Delta = 0$.

```{figure} /img/odes/perturbed_roots.png
:name: perturbed_roots
:alt: perturbed_roots
:align: center

Perturbation of repeated root (left) and special choice of coefficients (right)
```

**Second Method: Take the derivative with respect to $r$**

Note that the limit {eq}`ref-odes-review-second-12` is equivalent to the derivative

$$
\lim_{\epsilon \to 0} y(x,\epsilon) = \left. \frac{\partial}{\partial r} \left( e^{rx} \right) \ \right|_{r=r_1}
$$(ref-odes-review-second-14)

therefore we see that the new solution $x e^{r_1 x}$ was obtained by taking the derivative of $y(r,x)= e^{rx}$ with respect to $r$ and then making the substitution $r=r_1$. This is, in fact, a general procedure that we will use later in the course. To see why this procedure works, let

$$
y(r,x) = e^{rx}
$$(ref-odes-review-second-15)

and compute

$$
L \left[ y(r,x) \right] &= a(r - r_1)^2 e^{rx} \\
L \left. \left[ \frac{\partial y}{\partial r}(r,x) \right] \right|_{r=r_1} &= \left. \left( 2a(r - r_1) e^{rx} + a(r - r_1)^2 x e^{rx} \right) \right|_{r = r_1} = 0
$$(ref-odes-review-second-16)

Therefore

$$
\left. \frac{\partial y}{\partial r}(r,x) \right|_{r = r_1} = x e^{r_1 x}
$$(ref-odes-review-second-17)

is also a solution. Thus, to summarize, the general solution for the case of a repeated root is

$$
y(x) = c_1 e^{r_1 x} + c_2 x e^{r_1 x} \ \ , \ \ \ c_1,c_2 \in \mathbb{R}
$$(ref-odes-review-second-18)

## Case III: Complex Conjugate Roots

If $\Delta < 0$ then $r_1$ and $r_2$ are **complex conjugate roots**: $r_1,r_2 \in \mathbb{C}$ and $r_1 = \bar{r}_2$. Let

$$
\lambda = -\frac{b}{2a} \hspace{10mm} \text{and} \hspace{10mm} \mu = \frac{\sqrt{4ac-b^2}}{2a}
$$

such that 

$$
r_1 = \lambda + i \mu \hspace{10mm} \text{and} \hspace{10mm} r_2 = \lambda - i \mu
$$(ref-odes-review-second-19)

The general solution of {eq}`ref-odes-review-second-0` in this case is

$$
y(x) &= c_1 e^{(\lambda + i\mu)x} + c_2 e^{(\lambda - i\mu)x} \\
&= e^{\lambda x} \left( A \cos \mu x + B\sin \mu x \right) \ \ , \ \ \ A,B \in \mathbb{R}
$$(ref-odes-review-second-20)

## Examples

**Exmaple 1.** Find the general solution of the equation

$$
L[y] = y'' + y' - 6y = 0
$$

The roots of the characteristic polynomial are $r_1 = 2$ and $r_2 = -3$ and the general solution is

$$
y(x) = c_1 e^{-3x} + c_2 e^{2x} \ \ , \ \ \ c_1,c_2 \in \mathbb{R}
$$

**Exmaple 2.** Find the general solution of the equation

$$
L[y] = y'' + 6y' + 9y = 0
$$

The roots of the characteristic polynomial are repreated $r_1 = r_2 = -3$ and the general solution is

$$
y(x) = c_1 e^{-3x} + c_2 x e^{-3x} \ \ , \ \ \ c_1,c_2 \in \mathbb{R}
$$

**Exmaple 3.** Find the general solution of the equation

$$
L[y] = y'' - 4y' + 13y = 0
$$

The roots of the characteristic polynomial are $r_1 = 2 + 3i$ and $r_2 = 2 - 3i$ and the general solution is

$$
y(x) = e^{2 x} \left( A \cos 3x + B \sin 3x \right)
$$