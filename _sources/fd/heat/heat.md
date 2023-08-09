# FTCS Method

## Forward Time Central Space (FTCS)

Consider the following initial-boundary value problem for the heat equation

$$
\frac{\partial u}{\partial t} &= \alpha^2\frac{\partial ^2u}{\partial x^2} \quad 0<x<1, \quad t>0 \\
\text{BC:}\quad u(0,t) &= u(1,t) = 0 \\
\text{IC:}\quad u(x,0) &= f(x)
$$(HeatEq)

The basic idea is to replace the derivatives in the heat equation by difference quotients. We consider the relationships between $u$ at $(x,t)$ and its neighbours a distance $\Delta x$ apart and at a time $\Delta t$ later.

Corresponding to the difference quotient approximations introduced in {ref}`fd-intro`, we consider the following partial difference approximations.

Start with the Taylor series with respect to $t$

$$
u(x,t+\Delta t) = u(x,t)+\Delta t\frac{\partial u}{\partial t} (x,t)+\frac{\Delta t^2}{2!}\frac{\partial ^2u}{\partial t^2} (x,t)+\cdots
$$(ref-fd-pdes-9)

and then re-arrange and divide by $\Delta t$ to obtain the **forward difference in time** approximation

$$
\frac{u(x,t+\Delta t)-u(x,t)}{\Delta t} =\frac{\partial u}{\partial t}(u,t)+O(\Delta t)
$$(TimeDifference)

Now consider the following Taylor series with respect to $x$

$$
u(x + \Delta x,t) &= u(x,t) + \Delta x\frac{\partial u}{\partial x}(x,t) +\frac{\Delta x^2}{2!}\frac{\partial^2u}{\partial x^2}(u,t) \\
& + \frac{\Delta x^3}{3!}\frac{\partial^3u}{\partial x^3}(x,t) + \frac{\Delta x^4}{4!}\frac{\partial ^4u}{\partial x^2}(x,t)+\cdots \\
u(x - \Delta x,t) &= u(x,t)-\Delta x\frac{\partial u}{\partial x}(x,t) +\frac{\Delta x^2}{2!}\frac{\partial ^2u}{\partial x^2}(x,t) \\
& - \frac{\Delta x^3}{3!}\frac{\partial^3u}{\partial x^3}(x,t) +\frac{\Delta x^4}{4!}\frac{\partial ^4u}{\partial x^4}(x,t)+\cdots
$$(ref-fd-pdes-11)

Adding and re-arranging yields the **central differences in space** approximation

$$
\frac{u(x+\Delta x,t)-2u(x,t)+u(x-\Delta x,t)}{\Delta x^2}=\frac{\partial ^2u}{\partial x^2}(x,t) +O(\Delta x^2)
$$(CentralDiffSpace)

Substituting {eq}`TimeDifference` and {eq}`CentralDiffSpace` into {eq}`HeatEq` we obtain

$$
\frac{u(x,t+\Delta t)-u(x,t)}{\Delta t}=\alpha ^2\left( \frac{u(x+\Delta x,t)-2u(x,t)+u(x-\Delta x,t)}{\Delta x^2}\right) +O(\Delta t,\Delta x^2)
$$(ref-fd-pdes-13)

Re-arranging we get

$$
u(x,t+\Delta t)=u(x,t)+\alpha ^2\left(\frac{\Delta t}{\Delta x^2}\right) \left\{ u(x+\Delta x,t)-2u(x,t)+u(x-\Delta x,t)\right\}
$$(ref-fd-pdes-14)

We subdivide the spatial interval $[0,1]$ into $N+1$ equally spaced sample points $x_n=n\Delta x$. The time interval $[0,T]$ is subdivided into $M+1$ equal time levels $t_k=k\Delta t$. At each of these space-time sample points we introduce approximations:

$$
u(x_n,t_k) \simeq u_n^k.
$$(ref-fd-pdes-15)

```{figure} /img/fd/heat_eq_scheme.png
:name: fd_scheme
:alt: fd_scheme
:align: center
```

$$
u_n^{k+1}=u_n^k+\alpha ^2\left(\frac{\Delta t}{\Delta x^2}\right) \left(u_{n+1}^k-2u_n^k+u_{n-1}^k\right)
$$(ref-fd-pdes-16)

## Implementing Derivative Boundary Conditions

Assume that the boundary conditions for {eq}`HeatEq` are changed to

$$
\text{BC:}\quad u(0,t)=0,\quad\frac{\partial u}{\partial x}(1,t)=0
$$(ref-fd-pdes-17)

Consider a central difference approximation to $\displaystyle\frac{\partial u}{\partial x}(1,t)$, where $x_N=N \Delta x = 1$,

$$
\frac{u(x_N+\Delta x,t)-u(x_N-\Delta x,t)}{\Delta x}=0
$$(ref-fd-pdes-18)

Re-arranging we obtain:

$$
u(x_N+\Delta x,t)=u(x_N-\Delta x,t)
$$(ref-fd-pdes-19)

Since $x_N=1$ we observe that $x_N+\Delta x$ is outside the domain. To accomodate this we introduce an extra column $u_{N+1}$ into which we copy the values $u_{N-1}$. In the column $x_N$ we implement the same difference approximation for the Heat Equation, namely:

$$
u_N^{k+1}=u_N^k+\alpha ^2( \frac{\Delta t}{\Delta x^2} ) (u_{N+1}^k-2u_N^k+u_{N-1}^k)
$$(ref-fd-pdes-20)

While $u_{N+1}^k=u_{N-1}^k $ (see {ref}`ref-fd-pdes-19`) since column $u_{N-1}^k$ is copied to column $u_{N+1}^k$. 

Note that this BC could be implemented another way without introducing the additional column, by eliminating $u_{N+1}$ from $(*)$ and $(**)$:

$$
u_{N}^{k+1}=u_N^k+2\alpha ^2\left(\frac{\Delta t}{\Delta x^2}\right)\left( u_{N-1}^k-u_N^k\right)
$$(ref-fd-pdes-21)

If this latter equation is implemented at $x_N$ there is no need to introduce an extra column $u_{N+1}$ or to implement the difference equation given in (**) as the the derivative boundary condition is taken care of automatically.

## Stability of FTCS  

Consider the following finite difference approximation to the 1D heat equation:

$$
u_n^{k+1}-u_n^k=\frac{\Delta t}{\Delta x^2} \left( u_{n+1}^k-2u_n^k+u_{n-1}^k\right)\quad \mbox{where} \quad u_n^k\simeq u(x_n,t_k)
$$(ref-fd-pdes-22)

Let $\displaystyle u_n^k=\phi_ke^{in\Delta x\theta}$ then

$$
(\phi_{k+1}-\phi_k) e^{in\Delta x\theta} &= \frac{\Delta t}{\Delta x^2} \left( e^{i\Delta x\theta}-2+e^{-i\Delta x\theta}\right) \phi_k e^{in\Delta x\theta}\\
&= \frac{\Delta t}{\Delta x^2}\left[ 2\cos (\theta\Delta x)-2\right]\phi_k e^{in\Delta x\theta}
$$(ref-fd-pdes-23)

Therefore

$$
\phi_{k+1} &= \phi_k-\frac{\Delta t}{\Delta x^2}4\sin^2\left(\frac{\theta\Delta x}{2}\right)\phi_k \quad \mbox{since} \quad\cos (\theta\Delta x)-1=-2\sin^2 \left(\frac{\theta\Delta x}{2}\right)\\
&= \left[1-\frac{4\Delta t}{\Delta x^2}\sin^2\left(\frac{\theta\Delta x}{2}\right)\right]\phi_k
$$(ref-fd-pdes-24)

Now for stability we require that $|\phi_{k+1}|\leq |\phi_k|$ so that

$$
\left| 1-\frac{4\Delta t}{\Delta x^2}\sin^2\left(\frac{\theta\Delta x}{2}\right)\right|\leq 1
$$

and therefore

$$
-2\leq -\frac{4\Delta t}{\Delta x^2}\sin^2\left(\frac{\theta\Delta x}{2}\right)\leq 0
$$(ref-fd-pdes-25)

The right inequality is satisfied automatically, while the left inequality can be re-written  in the form:

$$
\frac{4\Delta t}{\Delta x^2}\sin^2\left(\frac{\theta\Delta x}{2}\right)\leq 2
$$(ref-fd-pdes-26)

Since $\sin^2(\theta \Delta x / 2) \leq 1$ for all $\theta$, the condition {ref}`ref-fd-pdes-26` is satisfied provided

$$
\Delta t\le\frac{\Delta x^2}{2}
$$(eqStabilityCondHeat)

## Exercises

**Numerical Instability:**

(a) Change the $\Delta t$ in cell D1 from $0{.}001$ to $0{.}05$ and you will observe what is known as a numerical instability. Now change $\Delta t$ to $0{.}00625$, which is known as the stability boundary predicted by {eq}`eqStabilityCondHeat` and observe what happens. Now let $\Delta t=0{.}006$ and observe the abrupt change in the solution - it is much closer to what we would expect.

(b) Derive the stability condition for the finite difference approximation of the 1D heat equation when $\alpha^2\ne 1$.

$$
u_n^{k+1}-u_n^k=\frac{\alpha^2\Delta t}{\Delta x^2}\left( u_{n+1}^k-2u_n^k+u_{n-1}^k\right)
$$(ref-fd-pdes-28)

**Truncation Error:**

The instability noted in 1 above is not the only source of error in the numerical approximation. Although numerical instability is evident for a parameter choice that is unstable, the other type of error is present in almost every type of numerical approximation scheme.  This class of error results from discarding the $O(\Delta x^2)$ and $O(\Delta t)$ terms in {eq}`CentralDiff1` and {eq}`ForwardDiff` when we replace derivatives in {eq}`HeatEq` by difference quotients. This error is known as the truncation error. To estimate the magnitude of the truncation error, change the spread sheet to implement the initial condition

$$
f(x)=\left\{\begin{array}{ll}2x &0<x<1/2\\ 2(1-x)&1/2\leq x <1\end{array}\right. .
$$(ref-fd-pdes-29)

Now code up the Fourier Series (in another spread sheet) that is derived in Lecture 10, Exercise 10.1 and compare the numerical solution to the "exact" Fourier Series solution with 50 terms. The difference between the two is mainly due to the truncation error since the round-off error is about $10^{-12}$ and does not grow if stable parameters are used.

**Derivative Boundary Conditions:**

Implement a derivative boundary condition the left endpoint $x=0$. Check the numerical solution against the problem solved in Lecture 11.