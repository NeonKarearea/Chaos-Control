### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 52716ab2-d273-11ed-1f86-09f8d71d4725
md"""
# Mid-term Report: CJ Barnes-Wilkie"""

# ╔═╡ 6994f507-6138-418a-bb79-f5338f1d8e0b
md"""
## Background and History"""

# ╔═╡ 8b1de72b-3f09-42e7-aeb7-0b6a3d5a8091
md"""
## Theory
For a single pendulum, we have a mass $m$ and a length $l$. We can release it from an angle $\theta$ and from there we can the position and velocity of the pendulum. We get:

$x = l\sin(\theta)$
$y = -l\cos(\theta)$
$\dot{x} = l\dot{\theta}\cos(\theta)$
$\dot{y} = l\dot{\theta}\sin(\theta)$

where $x$ is the $x$ position, $y$ is the $y$ position, $l$ is the length of the pendulum, $\theta$ is the angle, and the dotted variable are the time derivatives. For a single pendulum, it is easy to predict how it will behave, but what if we have two pendulums with one of them attached to the mass of the other? This particular system is much more difficult to solve by hand, as we can see by the position and velocity equations for the double pendulum. Below are the position equations.

$x_{1} = l_{1}\sin(\theta_{1})$
$y_{1} = -l_{1}\cos(\theta_{1})$
$x_{2} = x_{1} + l_{2}sin(\theta_{2}) = l_{1}\sin(\theta_{1})+l_{2}\sin(\theta_{2})$
$y_{2} = y_{1} - l_{2}cos(\theta_{2}) = -l_{1}\cos(\theta_{1})-l_{2}\cos(\theta_{2})$where variables with the 1 subscript are for the first pendulum, and variables with the 2 subscript are for the second variable. Below are the velocity equations.

$\dot{x_{1}} = l_{1}\dot{\theta_{1}}\cos(\theta_{1})$
$\dot{y_{1}} = l_{1}\dot{\theta_{1}}\sin(\theta_{1})$
$\dot{x_{2}} = \dot{x_{1}} + l_{2}\dot{\theta_{2}}\cos(\theta_{2}) = l_{1}\dot{\theta_{1}}\cos(\theta_{1}) + l_{2}\dot{\theta_{2}}\cos(\theta_{2})$
$\dot{y_{2}} = \dot{y_{1}} + l_{2}\dot{\theta_{2}}\sin(\theta_{2}) = l_{1}\dot{\theta_{1}}\sin(\theta_{1}) + l_{2}\dot{\theta_{2}}\sin(\theta_{2})$

From here we can follow the methodology of the paper *Chaos from Simplicity: An Introduction to the Double Pendulum* by Chen, J. to get the necessary equations of motion. Those equations have been put into the **Numerical Methods** section.

Another thing we can look at for chaotic systems is the Lyapunov exponential. We can  find this from the eventual $\omega_{(t)}$, so:

$\delta_0 = \omega_{(t_0 + \delta t)}-\omega_{(t_0)}$From this we can write:

$|\delta| = |\delta_0|e^{n\sigma}$where $\sigma$ is equal to the Lyapunov exponential. The Lyapunov exponential can be determined by adapting the formula from the Chaos-II lecture notes.

$\sigma = \lim_{n\to\infty}\frac{1}{n}\left\{\sum_{i=0}^{n-1}\ln\big|\dot{\omega}_{(t_i)}\big|\right\}$This allows us to determine if there is chaotic motion. If $\sigma$ is negative, then the system is regular, but if the $\sigma$ is positive, then the system is chaotic.
"""

# ╔═╡ bf90d8ba-b78e-4882-8325-277060ab96e7
md"""
## Numerical Methods
From *Chen, J.* we have the following equations:

$\omega_{1} = \dot{\theta_{1}}$

$\omega_{2} = \dot{\theta_{2}}$

$\dot{\omega_{1}} = \frac{m_{2}l_{2}\omega_{1}^{2}\sin(2\Delta\theta)+2m_{2}l_{2}\omega_{2}^{2}\sin(\Delta\theta)+2gm_{2}\cos(\theta_{2})\sin(\Delta\theta)+2gm_{1}\sin(\theta_{1})}{-2l_{1}(m_{1}+m_{2}\sin^{2}(\Delta\theta)}$

$\dot{\omega_{2}} = \frac{m_{2}l_{2}\omega_{2}^{2}\sin(2\Delta\theta)+2(m_{1}+m_{2})l_{1}\omega_{1}^{2}+2g(m_{1}+m_{2})\cos(\theta_{1})\sin(\Delta\theta)}{2l_{2}(m_{1}+m_{2}\sin^{2}(\Delta\theta))}$

where $m_{1}$ is the mass of the first pendulum's bob, $m_{2}$ is the mass of the second pendulum's bob, and $g$ is the acceleration due to gravity ($9.81\;m\;s^{-1}$). From here we can set up the function that contains all of the equations. From here, we can use the package DifferentialEquations.jl and solve over a particular time range. We will be testing 4 algorithms: Vern6, Vern6 with the lazy interpolant off, Vern7, and Vern7 with the lazy interpolant off. This means that we can figure out which one is more accurate, and it also allows us to compare with the algorithm that Chen used (The method used was ode113 in MATLAB, which equates to the VCABM() algorithms, which is a less efficient version of Vern7). Then we can use ProblemODE to set up the problem and the inital conditions for the computer to solve, and then we can use the solve function to solve it for us.

There are two main checks we can run to see if out model works, those are the following:

	- Conservation of Energy: As this is an ideal system, the energy should be conserved. This means that if we calculate the kinetic and potential energy as each time step we expect the energy to be consistent. With this we can also test the 4 different algorithms and find the one that conserves the energy the best.

	- Similarity to a simple pendulum: By setting both of the angles small, we expect the double pendulum to act like a regular pendulum.

Once we have that, we need to determine the Lyapunov exponential. We can use the package ChaosTools.jl to help with that. We can use the function ChaosTools.lyapunovspectrum to find the range of possible lyapunov exponentials. We can also just use ChaosTools.lyapunov to get the largest exponent. We can also see if the Lyapunov exponents are correct by checking where the regular section is. This should fit the simple pendulum.
"""

# ╔═╡ Cell order:
# ╟─52716ab2-d273-11ed-1f86-09f8d71d4725
# ╟─6994f507-6138-418a-bb79-f5338f1d8e0b
# ╟─8b1de72b-3f09-42e7-aeb7-0b6a3d5a8091
# ╟─bf90d8ba-b78e-4882-8325-277060ab96e7
