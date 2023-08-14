---
title: First-Order Partial Differential Equations
author: Daning Huang
date: 07/16/2023
header-includes:
  - \usepackage{amsmath}
  - \usepackage{xcolor}
---

\newcommand{\dd}{\mathrm{d}}
\newcommand\ppf[2]{\frac{\partial #1}{\partial #2}}
\newcommand\ddf[2]{\frac{\dd #1}{\dd #2}}
\newcommand\Cr[1]{{\color{red} #1}}
\newcommand\Cg[1]{{\color{green} #1}}
\newcommand\Cb[1]{{\color{blue} #1}}

<mark>*Note: This chapter is expected to be after two chapters on second-order PDEs.  So the students are expected to know some basics about PDEs, including a linear version of method of characteristics.*</mark>

<mark>Heroku in the text refers to examples on my website https://aersp313.herokuapp.com/, as we have not converged on how to integrate the interactions.</mark>

# Introduction
In this chapter we consider a special type of Partial Differential Equations (PDEs) that have only first-order partial derivatives, i.e., first-order PDEs.  Just like the second-order PDEs that we have seen, first-order PDEs play a fundamental role in many engineering applications.

One important example is the Euler equation for compressible inviscid flow; aerodynamic solvers based on Euler equations have been widely used as cost-effective tools for the preliminary design of subsonic, transonic and supersonic aircraft.  The simplified one-dimensional Euler equation reads
$$
\begin{aligned}
\rho_t + u\rho_x &= 0 \\
u_t + \frac{1}{\rho}p_x &= 0 \\
e_t + ue_t &= 0
\end{aligned}
$$
where $\rho$ is the gas density, $u$ the velocity, and $e$ the energy; $p$ is pressure as a function of $\rho$ and $u$.  As usual, the subscripts denote partial derivatives, e.g., $u_t = \ppf{u}{t}$.  Clearly, the Euler equations are PDEs that only involve first-order derivatives.

Another example is the traffic flow equation,
$$
u_t+uu_x=0,\quad u(x,0)=\phi(x)
$$
where $u$ is the "density" of vehicles, $\phi$ gives the initial distribution of vehicles, and the PDE predicts how the vehicles would be distributed over time.  Unlike the previous example, the traffic flow equation is a scalar one, i.e., there is only one unknown variable.  Yet, this surprisingly simple equation can explain the so-called phantom traffic jams, where dense traffic comes to a halt and restarts for no apparent reason.

# Model Problem
This chapter will focus on the scalar first-order PDEs, and solve a model problem in the following form:
$$
\begin{aligned}
\text{PDE:} &\quad u_t+V(u,x,t)u_x = S(u,x,t) \\
\text{Initial condition:} &\quad u(x,0) = \phi(x)
\end{aligned}
$$ {#eq:model}
where $u(x,t)$ is defined for $x\in(-\infty,\infty)$ and $t\in[0,\infty)$, i.e., the entire upper part of $x-t$ plane.  Due to the presence of the initial condition, that are defined on the entire $x$-axis, [@Eq:model] is an **initial value problem** (IVP).

## What does the PDE really mean? {#sec:mean}
We got three terms in the equation and let's look at them one by one.  To facilitate a qualitative and tangible understanding, we can think $u$ as the density distribution of some gas that varies over space $x$ and time $t$; in this setup, the initial density distribution is given by $\phi(x)$.  Also, for simplicity, assume $V$ and $S$ are constant.

First, the $u_t$ term is perhaps the most straightforward one to understand - it means the rate of change in the unknown variable, e.g., how the density vary over time.  Next, ignoring $V u_x$, let's look at the combined effect of $u_t$ and $S$, i.e., $u_t=S$; this equation simply defines the rate of change in $u$.  In our example, $S$ functions as a source (sink) that adds (removes) gas from the density distribution.  Hence, we usually call $s$ the **source** term, with the understanding that a sink is as if a "negative source".

Subsequently, let's look at the combined effect of $u_t$ and $V u_x$, ignoring the effect of the source term $S$.  This gives us an equation
$$
u_t + V u_x = 0,\quad u(x,0) = \phi(x)
$$ {#eq:meaning}
When $V$ is constant, we can verify that the solution to [@Eq:meaning] is
$$
u(x,t) = \phi(x-Vt)
$$ {#eq:msol}




> Exercise: Verify the solution as an exercise.  This can be done by checking if [@Eq:msol] satisfies *both* the PDE and initial condition in [@Eq:meaning].
> 
> We will show how to find this solution very soon.

What does [@Eq:msol] mean?  If we substitute in some numbers: At $t=1$, $u(x,t) = \phi(x-V)$, and $\phi(x-V)$ is simply the initial condition $\phi(x)$ shifted to the right by $V$; then at $t=2$, $u(x,t) = \phi(x-2V)$ is $\phi(x)$ shifted to the right by $2V$, and so on.  To summarize, [@Eq:msol] is just representing that the initial condition $\phi(x)$ is constantly shifted to the right at a rate of $V$, see **Visualization**.  In our example, this can mean that there is constantly a wind blowing from left to right at the velocity $V$, that moves the initial gas distribution to the right at the same rate; this phenomenon is **advection**.

At this point, we can conclude that, qualitatively, the first-order PDE describes the time variation of some quantity due to the combined effects of *source* and *advection*.

<mark>Heroku First-order PDE (3D) Case 1</mark>

## Plan of attack
Now return to the more general model problem [@Eq:model].  You might have noticed two key differences from the second-order PDEs we have solved so far:

1. [@Eq:model] does not have boundary conditions.  There are fewer conditions to satisfy and our life is thus *easier*.
2. [@Eq:model] can be nonlinear, i.e., we might have terms involving products of unknown functions and their derivatives.  For example, if $V(u,x,t)=u^2$, we would get a nonlinear convection term $u^2 u_x$ in the PDE.  The nonlinearity would make our life *harder*.

Yet, luckily, [@Eq:model] can be viewed as *quasi-linear*, in the sense that the $u_t$ and $u_x$ terms are always first-order regardless of the nonlinearity in $V$ and $S$.  Here "first-order" means that we do not have "strange" terms such as $u_t^2$, $\sqrt{u_x}$, $\sin(u_t)$, etc.  We will leverage this quasi-linearity and use a modified version of the **method of characteristics** to solve the IVP of first-order PDEs.

# Method of Characteristics: Nonlinear Version
## What are characteristics

### Definition
Let's look at again the simple problem of advection with constant velocity in [@Sec:mean]:
$$u_t+Vu_x=0,\quad u(x,0)=\phi(x)$$
the solution is [@Eq:msol], $u(x,t)=\phi(x-Vt)$.  If we define a new variable $\eta(x,t)=x-Vt$, then $u$ can be written as a *uni-variate* function $u(\eta)$.  This fact actually has a non-trivial implication for constant advection: For whatever combination of $(x,t)$, as long as $\eta$ is the same, $u$ remains the same.

<mark>Heroku First-order PDE (2D) Case 1 Plot 3</mark>

<iframe src="myplot_small.html" width="600" height="400">
</iframe>


Let's use the **Visualization** to get a geometrical intuition of what the implication means.  In the numerical example, we simply set $V=1$, and $\eta=x-Vt=x-t$ represents a straight line with a slope $1/V=1$ that intersect $x$-axis at $x=\eta$.   Also, we pick an initial condition for the IVP
$$
\phi(x)=\exp(-(x-2)^2/2)
$$
This function is like a single wave: It has a peak of 1 at $x=2$ and is almost zero when $x<1$ or $x>3$.

Now we start from an arbitrary point on the $x$-axis, say $x=\eta=2$; at this point $\phi(2)=1$ and $t=0$.  Associated with the point $(x,t)=(2,0)$, we have a curve $2=x-t$.  We can verify that the values of $u$ at all the points on this curve, e.g., $(3,1)$, $(4,2)$, $(5,3)$, etc., are the same as the value at the point on the $x$-axis $(x,t)=(2,0)$, i.e.,
$$
u(3,1)=u(4,2)=u(5,3)=\cdots=u(2,0)=\phi(2)=1
$$
Furthermore, note that the choice of point on the $x$-axis is arbitrary.  This means we can just pick any $x=\eta$ on the $x$-axis, and find a curve $\eta=x-Vt$; on this curve $u$ is constantly $\phi(\eta)$.

> Exercise: Consider $\phi(1.5)\approx 0.5$, and identify the combinations $(x,t)$ for which $u=\phi(1.5)$.

To sum up at this point, in the constant advection problem, the equation $\eta=x-Vt$ represents a family of special curve with a parameter $\eta$; every point on the $x$-axis emanates one curve and it runs forward in time.  On each curve the variation of $u$ is *greatly simplified*, in the current case of which $u$ becomes constant.  This type of special curve is called the characteristic lines, or **characteristics** for short.

### Fixed and moving frames
So far we have been looking at one characteristics at a time; this is as if we slice the solution $u$ and the $(x,t)$ plane along the characteristics, so that we only look at this slice.

There are certainly other ways to slice the solution $u$, one of which is to slice at a time $t$.  In this case, we can think the characteristics as *rails*.  Along the rails, a slice of solution at $t=0$, $u(x,0)$, moves forward by time $t$, and arrive at another slice of solution $u(x,t)$.  Depending on the shape of rails, the slice can shift to right or left, dilute or shrink, and so on.  In the constant advection problem, the characteristics run straight to the right, hence the slice of solution simply shift to the right at constant rate.  We call this perspective of "characteristics as rails" as **fixed frame**, since we are as if sitting on the $x$-axis and watch the slice of solution moving forward in time.

There is also a view of **moving frame**.  Following the analogy of rails, then we ourselves are simply moving on the rails together with the solution.  As a result, we would not sense any shift in the slice of solution, and rather would sense only the change in the distribution of the solution.  In the constant advection problem, the shape remains the same, hence in the moving frame, we would sense *no* change at all.  Mathematically, the moving frame is equivalent to have a change of variables and write the solution in terms of $\eta$ and $t$.

<mark>Heroku First-order PDE (2D) Case 1 Plots 1 and 2</mark>

## Finding characteristics
The key takeaway from the previous section is, if we can find the characteristics of a first-order PDE, then we probably can use the characteristics to simplify and solve the IVP.  Now the question becomes how to find the characteristics.  To do so we will leverage the concept of total derivative from multi-variate calculus.

The total derivative of a multivariate function $u(x,t)$ with respect to $t$ is,
$$
\Cr{\ddf{u}{t}} = u_x \Cb{\ddf{x}{t}} + u_t
$$ {#eq:td}
Now we rearrange the first-order PDE and compare it side by side with [@Eq:td],
$$
\Cr{S} = u_x \Cb{V} + u_t
$$ {#eq:pde}
Clearly, there is a striking similarity between [@Eq:td] and [@Eq:pde], highlighted by the color coding.  Since the total derivative formula [@Eq:td] is valid for an arbitrary function, the similarity indicates that our first-order PDE should be representing some form of total derivative of $u$ *as well*.  In other words, it seems natural to set
$$
\ddf{u}{t} = S,\quad \ddf{x}{t}=V
$$ {#eq:odes}
If we solve the second ODE starting from a point $P$ on the $x$-axis, e.g., $x=\eta$, $t=0$, then we would find $x$ as a function of $t$ that passes through the point $(\eta,0)$.  For example, in the constant advection problem, and we can solve
$$
\ddf{x}{t}=V,\quad x(t=0)=\eta
$$
to find $x=Vt+\eta$, or $\eta=f(x,t)=x-Vt$.  This is exactly what we called as characteristics, and this line starts from the point $(\eta,0)$ and runs forward in time.

Next, at point $P$, we also know $u(\eta,0)=\phi(\eta)$, so we can supply the first ODE in [@Eq:odes] with initial conditions $t=0$, $u=\phi(\eta)$.  For example, in the constant advection problem, we have $S=0$ and
$$
\ddf{u}{t}=0,\quad u(t=0)=\phi(\eta)
$$
Clearly the solution is $u=\text{const}=\phi(\eta)$, which is valid along the characteristics.

Lastly, once we have found the characteristics equation $\eta=f(x,t)$ and the solution $u$ along the characteristics, we can combine the solutions to find the solution to the first-order PDE.  Again, for the example of constant advection, we know $\eta=x-Vt$ and $u=\phi(\eta)$, and thus the final solution is
$$
u(x,t) = \phi(\eta) = \phi(x-Vt)
$$
just exactly as what we found earlier in [@Eq:msol].

## Summary of method
The above discussion lays out the steps of the method of characteristics for the IVP of a first-order PDE.  In the following we summarize the method as a **trilogy**

1. Convert the problem into two IVPs of ODEs:
$$
\begin{aligned}
\ddf{x}{t}&=V,\quad x(t=0)=\eta \\
\ddf{u}{t}&=S,\quad u(t=0)=\phi(\eta)
\end{aligned}
$$
1. Solve the two ODEs and obtain (1) $x$ as a function of $(t,\eta)$ and (2) $u$ as a function of $(t,\eta)$.
2. Rearrange $x(t,\eta)$ to obtain the characteristics $\eta=f(x,t)$; eliminate $\eta$ in $u$ and obtain the final solution $u(x,t)$.

# Examples
Here we provide a series of examples of increasing complexity, and illustrate the typical behaviors of first-order PDEs.  When not explicitly specified, we continue to use the single wave equation, $\phi(x)=\exp(-(x-2)^2/2)$, as the initial condition for the IVP.

## Non-uniform advection: Homogeneous case
We consider the following problem
$$
u_t+e^{-t}u_x = 0,\quad u(x,0) = \phi(x)
$$ {#eq:eg1}

First, we identify that $V=e^{-t}$ and $S=0$, which gives us the ODEs
$$
\begin{aligned}
\ddf{x}{t}&=V=e^{-t},\quad x(t=0)=\eta \\
\ddf{u}{t}&=S=0,\quad u(t=0)=\phi(\eta)
\end{aligned}
$$
Here $\eta$ means we are picking an arbitrary point $x=\eta$ on the $x$-axis.

Next, we solve the first ODE for the characteristics.  From the ODE we find $x = -e^{-t} + d_0$ with undetermined coefficient $d_0$.  From the IC, $\eta=-e^{-0}+d_0$, we find $d_0=\eta+1$, so the characteristics is
$$
x = -e^{-t} + \eta + 1,\quad\text{or}\quad \eta = x + e^{-t} - 1
$$
Intuitively, thinking $\ddf{x}{t}$ as the speed of advection, the ODE indicates that the speed decreases over time.  We can confirm the intuition by plotting the characteristics for several values of $\eta$.  Indeed, initially the characteristics move to the right but gradually come to a halt as time increases.

Subsequently, we solve the second ODE that tells us the evolution of $u$ along one characteristics.  In this example, the ODE is easy to solve, $u$ stays constant all the time:
$$
u=\phi(\eta)
$$
Intuitively, since the source term is zero, nothing is added or removed from $u$, and hence $u$ is not changing over time.

Lastly, we combine the ODE solutions to obtain the complete PDE solution:
$$
u(x,t)=\phi(\eta)=\phi(x + e^{-t} - 1)
$$
Looking at the visualization, the characteristics are as if rails and the initial condition moves exactly along the characteristics.  And again, imagine if we are in the moving frame of the characteristics, i.e., moving at the given advection speed, then we would see $u$ stays absolutely unchanged starting from the initial condition.

<mark>Heroku First-order PDE (2D) Case 2</mark>

## Non-uniform advection: Non-homogeneous case
Next, we modify the problem a little bit
$$
u_t+e^{-t}u_x = -u^2,\quad u(x,0) = \phi(x)
$$ {#eq:eg2}

First, we identify that $V=e^{-t}$ and $S=-u^2$, which gives us the ODEs
$$
\begin{aligned}
\ddf{x}{t}&=V=e^{-t},\quad x(t=0)=\eta \\
\ddf{u}{t}&=S=-u^2,\quad u(t=0)=\phi(\eta)
\end{aligned}
$$

Next, the first ODE is the same as the previous example, and the characteristics has been found to be
$$
\eta = x + e^{-t} - 1
$$

Subsequently, we solve the second ODE, that is more complex than the previous example.  From the ODE,
$$
\begin{aligned}
\frac{\dd u}{-u^2} &= \dd t \\
\frac{1}{u} &= t + d_0 \\
u &= \frac{1}{t+d_0}
\end{aligned}
$$
From the IC,
$$
u(t=0) = \phi(\eta) = \frac{1}{0+d_0}\quad\Rightarrow\quad d_0=\frac{1}{\phi(\eta)}
$$
so the final ODE solution is
$$
u = \frac{1}{t+d_0} = \frac{1}{t+1/\phi(\eta)} = \frac{\phi(\eta)}{t\phi(\eta)+1}
$$

Lastly, we combine the ODE solutions to obtain the complete PDE solution:
$$
u(x,t)=\frac{\phi(x + e^{-t} - 1)}{t\phi(x + e^{-t} - 1)+1}
$$
Looking at the visualization, the characteristics are still functioning as if guidelines and the initial condition moves exactly along the characteristics.  But the difference is that $u$ itself is changing as it moves along the characteristics.  If we look at the moving frame, while $u$ is not shifting to right or left as usual, its amplitude decreases as time increases.

<mark>Heroku First-order PDE (2D) Case 3</mark>

We can also compare side by side the solutions of the homogeneous and non-homogeneous cases in 3D surface plots.  The two plots show the same curved distortion to the positive $x$-direction as time increases, which is due to the identical set of characteristics.  The difference is that, the amplitude of solution in homogeneous case remains constant over time, while the amplitude in the other case gradually decays; the difference is attributed to the extra source term in the latter case.

> Exercise: What if the source term is $u^2$?  First think intuitively, and then verify your intuition.

<mark>Heroku First-order PDE (3D) Cases 2 and 3, side by side</mark>

## Nonlinear advection
Lastly, we touch upon a more complex problem,
$$
u_t+uu_x = 0,\quad u(x,0) = \phi(x) = 1-x^2
$$ {#eq:eg3}
You might recognize this equation from the motivating example of traffic flow.

We again start with identifying the ODEs
$$
\begin{aligned}
\ddf{x}{t}&=u,\quad x(t=0)=\eta \\
\ddf{u}{t}&=0,\quad u(t=0)=1-\eta^2
\end{aligned}
$$
Note that now the ODE for characteristics basically says that the advection speed is larger when $u$ is larger, and the speed may become negative if $u$ is negative; the ODE now involves the unknown function $u$ and cannot be solved directly.

Instead, we attempt to first solve the ODE for $u$, which, luckily, is rather simple:
$$
u=\text{const}=1-\eta^2
$$
Subsequently, returning to the first ODE, which now becomes,
$$
\ddf{x}{t}=u=1-\eta^2,\quad x(t=0)=\eta
$$
and we have $x=(1-\eta^2)t + d_0$.  Using the IC, $\eta=(1-\eta^2)\times 0 + d_0$, we get $d_0=\eta$, so the characteristics equation is
$$
x=(1-\eta^2)t + \eta
$$
The first complication we see here is that $\eta$ is an implicit function of $(x,t)$.  But still, we could find an explicit expression of $\eta$ by solving a quadratic equation, and find
$$
\eta = \frac{1-\sqrt{1+4t^2-4xt}}{2t}
$$

> Exercise: A quadratic equation always has two roots.  Think about why we throw away one of the roots here.
>
> Hint: Look at $t=0$.

Lastly, we can combine the two ODE solutions to find the PDE solution,
$$
u(x,t) = 1-\eta^2 = \frac{-1+2tx+\sqrt{1+4t^2-4xt}}{2t^2}
$$ {#eq:eg3s}
The solution looks relatively more complex than previous ones.  The visualization shows the second complication.  It is clear that this time the characteristics are no longer parallel to each other.  In some regions, the lines move away from each other, while in some other regions the lines move closer.  This non-uniformity is due to the dependency of advection speed on the solution $u$ itself.  As a result, while $u$ would still move along the characteristics, the distribution would be distorted by the characteristics.

<mark>Heroku First-order PDE, a taste of shocks (2D)</mark>

> Problem solved?  Is [@Eq:eg3s] as "harmless" as it looks?  For example, is the solution defined when $t=1/2$ and $x>1$?
>
> This is where we see a "shock" in the solution.  We will discuss more in the advanced topics.

# Advanced Topics

<mark>*These will not be covered in class, and the full text will be developed in the final version of book.*</mark>

## More general case
<mark>*We give a more general treatment of the PDE.*</mark>
$$
\begin{aligned}
\text{PDE:} &\quad a(u,x,t)u_t+b(u,x,t)u_x = c(u,x,t) \\
\text{IC:} &\quad u(\xi) = u_0(\xi)\text{ on }x=x_0(\xi),\ y=y_0(\xi)
\end{aligned}
$$

## Existence of solution
<mark>*For mathematical rigor, we discuss when the solution exists.*</mark>

## When things break down
<mark>Now we formally formulate the problems of shocks and expansion waves, and discuss how to tackle these cases.</mark>
