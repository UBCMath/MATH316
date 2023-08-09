# Overview and Formulas 

## What is the Finite Difference Method?

*Under construction*

(fd-intro)=
## Finite Difference Formulas

Recall that the derivative of a function was defined by taking the limit of a difference quotient

$$
f'(x) = \lim_{\Delta x \to 0} \frac{f(x + \Delta x) - f(x)}{\Delta x}
$$(Derivative)

Now to use the computer to solve differential equations we go in the opposite direction - we replace derivatives by appropriate difference quotients. If we assume that the function can be differentiated many times then Taylor's Theorem is a very useful device to determine the appropriate difference quotient to use. For example, consider

$$
f(x + \Delta x) = f(x) + \Delta x f'(x) + \frac{\Delta x^2}{2!}f''(x) + \frac{\Delta x^3}{3!}f^{(3)}(x)+\frac{\Delta x^4}{4!}f^{(4)}(x) + \cdots
$$(Taylor)

Re-arranging terms in {eq}`Taylor` and dividing by $\Delta x$ we obtain

$$
\frac{f(x + \Delta x) - f(x)}{\Delta x} = f'(x) + \frac{\Delta x}{2}f''(x) + \frac{\Delta x^2}{3!}f^{(3)}(x) + \cdots
$$(ref-fd-pdes-2)

If we take the limit $\Delta x \to 0$ then we recover {eq}`Derivative`. But for our purposes it is more useful to retain the approximation

$$
\frac{f(x + \Delta x) - f(x)}{\Delta x} = f'(x) + \frac{\Delta x}{2}f''(\xi) = f'(x)+O(\Delta x)
$$(ForwardDiff)

We retain the term $\displaystyle\frac{\Delta x}{2}f''(\xi )$ in {eq}`ForwardDiff` as a measure of the error involved when we approximate $f'(x)$ by the difference quotient $\big( f(x+\Delta x)-f(x)\big) /\Delta x$. Notice that this error depends on how large $f''$ is in the interval $[x,x+\Delta x]$ (i.e. on the smoothness of $f$) and on the size of $\Delta x$. Since we like to focus on that part of the error we can control we say that the error term is of the order $\Delta x$ -- denoted by $O(\Delta x)$. Technically a term or function $E(\Delta x)$ is $O(\Delta x)$ if

$$
\frac{E(\Delta x)}{\Delta x} \stackrel{\Delta x\rightarrow 0}{\rightarrow} \quad \text{constant}
$$(ref-fd-pdes-4)

Now the difference quotient {eq}`ForwardDiff` is not the only one that can be used to approximate $f'(x)$. Indeed, if we consider the expansion of $f(x-\Delta x)$:

$$
f(x - \Delta x) = f(x) - \Delta xf'(x) + \frac{\Delta x^2}{2!}f''(x) - \frac{\Delta x^3}{3!}f^{(3)}(x) + \frac{\Delta x^4}{4!}f^{(4)}(x) + \cdots
$$(Taylorminus)

and we subtract {eq}`Taylorminus` from {eq}`Taylor`, and divide by $(2\Delta x)$ we obtain:

$$
\frac{f(x + \Delta x) - f(x - \Delta x)}{2\Delta x} = f'(x) + \frac{\Delta x^2}{3!}f^{(3)}(\xi)
$$(CentralDiff1)

We notice that the error term associated with this form of difference approximation is $O(\Delta x^2)$, which converges more rapidly to zero as $\Delta x\rightarrow 0$.

In order to obtain an approximation to $f''(x)$ we add {eq}`Taylorminus` to {eq}`Taylor`, which upon re-arrangement and dividing by $\Delta x^2$ leads to

$$
\frac{f(x + \Delta x) - 2f(x) + f(x - \Delta x)}{\Delta x^2} = f''(x) + \frac{1}{12}\Delta x^2f^{(4)}(\xi)
$$(CentralDiff2)

Due to the symmetry of the difference approximations {eq}`CentralDiff1` and {eq}`CentralDiff2` about the expansion point $x$ these are called central difference approximations. The difference approximation {eq}`ForwardDiff` is known as a forward difference approximation. We note that the central difference schemes {eq}`CentralDiff1` and {eq}`CentralDiff2` are second order accurate while the forward difference scheme {eq}`ForwardDiff` is only accurate to $O(\Delta x)$.