### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 34ac14f9-752a-4095-866f-0f32befc6957
using PlutoUI

# ╔═╡ 814f4a90-d274-11ed-346c-d7fcaa9dfadb
md"""
# Chaos Theory: James Elford
"""

# ╔═╡ ec7c6b22-1376-47d9-9b38-97ed2f420207
md"""
## History
"""

# ╔═╡ ba8bb4c6-be71-4adf-ba4f-c2f06dc431aa
md"""
Chaos theory is an area of mathmatics which began in the mid-1600's when Newton discoved differntial equations [1], in which we study systems which have chaotic properties. These properties include: a sensitivity to initial conditions, strange attractors, and non-periodicity. An early example of a system which involves chaos is the Three Body Problem, which is quite sensitive to both the initial positions and initial velocities of each body of mass. We should also note that there has been no general solution found for the Three Body Problem so far. [2]

One thing we can notice is that chaos can come from a very simple system, or a natural continuation of a simple system. We can see examples of this in problems such as the three body problem, which is a natural progression of the two body problem, the double pendulum which is an extension of the simple pendulum, and finally the logistic map which is a very simple looking equation. 
"""

# ╔═╡ 11825f7a-57de-4f9c-9000-f09a96736e44
md"""
## Backround
"""

# ╔═╡ f80d5bd6-fb20-4141-8ba8-99c2b6d1f790
md"""
As a starting point for our investigation into chaos, we will analyse a simple system, the logistic map, which is given by: 

$$x_{n+1} = rx_n(1 - x_n)$$

where $r \in [0,4]$.This appears like a very simple quadratic equation and is used to find the ratio between an existing population to some maximum population.

Firstly, we can easily find the fixed points by setting $x_n,x_{n+1} = x$, so $x = r - 1, x = 0$ are our fixed points

"""

# ╔═╡ 01788aa6-486f-4bee-941c-c42a4e2a5f8e


# ╔═╡ 0a7010e7-c093-4c26-afd3-cf4c17af7049
md"""
## Theory

"""

# ╔═╡ a0664415-813f-4e35-a519-ab44437720b9
double_pendulum_url = "https://www.researchgate.net/profile/Iman-Izadgoshasb/publication/331024972/figure/fig1/AS:725254309699586@1549925510035/Double-pendulum-system.png"

# ╔═╡ 46304beb-13f2-45d3-aef5-e9933472d42c
md"""$(Resource(double_pendulum_url))

Image taken from researchgate.net
"""

# ╔═╡ 3190dfb9-5bb3-41c9-81ae-1fcae71cd0a7
md"""
From an inspection of the double pendulm system, we can see that the equations of motion for the pendulums would depend on both the positions and velocities of each pendulum. From [2] we have the equations of motion as follows. 

!!! note "Positions"
	$x_1 = l_1 sin θ_1$
	$y_1 = -l_1cosθ_1$
	$x_2 = x_1 + l_2 sinθ_2 = l_1 sin θ_1 + l_2 sin θ_2$
	$y_2 = y_1 - l_2 sinθ_2 = -l_1 sin θ_1 - l_2 sin θ_2$
 
!!! note "Velocities"
	$ẋ_1 = l_1 \dot{θ_1} sin θ_1$
	$\dot{y_1} = l_1 \dot{θ_1} sin θ_1$
	$\dot{x_2} = \dot{x_1} + l_2 \dot{θ_2} sinθ_2 = l_1 \dot{θ_1} sin θ_1 + l_2\dot{θ} sin θ_2$
	$\dot{y_2} = \dot{y_1} + l_2 \dot{θ_2} sinθ_2 = l_1 \dot{θ_1} sin θ_1 + l_2\dot{θ} sin θ_2$

Where, $l_1,l_2$ are the lengths of each pendulum respectively $θ_1,θ_2$ are angles relative to the vertical. 

Then introducing the angular frequencys and letting $Δθ = θ_1 - θ_2$ reduces our system down to just 4 equations. 

!!! warning "Angular Frequencies"
	$ω_1 = \dot{θ_1}$
	$ω_2 = \dot{θ_2}$
	$\dot{ω_1} = \frac{m_2l_1ω_1^2sin(2Δθ) + 2m_2l_2ω_2^2sin(Δθ) + 2gm_2cosθ_2 sin(Δθ + 2gm_1sinθ_1)}{-2l_1(m_1+m_2sin^2Δθ)}$
	$\dot{ω_2} = \frac{m_2l_2ω_2^2sin(2Δθ) + 2(m_1+m_2)l_1ω_1^2sin(Δθ) + 2g(m_1 + m_2)cosθ_1sin(Δθ)}{2l_1(m_1+m_2sin^2Δθ)}$

Note that g is acceleration due to gravity. 

"""



# ╔═╡ b036d4ca-d8ea-43e1-ab23-bb99bbf91cec
md"""
## Numerical Methods
"""

# ╔═╡ f0c29b93-c63c-47cf-bcbd-9b7b58cb7e99
md"""
In order to analyse the chaos of the double pendulum we will be finding the Lyapunov 

We will be using the OrdinaryDiffEq, Plots, and chaotictools julia packages to help deal with the double pendulum system. 

We have also devised two methods in order to check our solutions are consistent with the laws of physics. As this is an idealised system, we know that energy will be conserved, so we can check the potential and kinetic energy at the initial conditions and at each time step. This will give us a way to quantify the amount of error accumulated by the ODE solving algorithm. 


Secondly, we can compare our solutions to that of a simple single pendulum when the angles are small, we know that it should behave like a single pendulum with length $l_1 + l_2$ 
"""

# ╔═╡ 9b2379cd-7dff-4f5d-9c05-0e7ef4a707a0


# ╔═╡ 55e85624-6196-4501-b04a-f259aeb84752
md"""
## Sources
[1] "Nonlinear Dynamics and Chaos With Applications to Physics, Biology, Chemistry,and Engineering" by Steven H. Strogatz

[2] 'The Three-Body Problem and the Equations of Dynamics' by Henri Poincaré, translated by Bruce D Popp. 

[3] 'Chaos from Simplicity: An Introduction to the Double Pendulum' by Joe Chen. 
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.50"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "d8b0bbb312600ec81f2769bd72048a77429debd9"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═814f4a90-d274-11ed-346c-d7fcaa9dfadb
# ╠═ec7c6b22-1376-47d9-9b38-97ed2f420207
# ╠═ba8bb4c6-be71-4adf-ba4f-c2f06dc431aa
# ╠═11825f7a-57de-4f9c-9000-f09a96736e44
# ╠═f80d5bd6-fb20-4141-8ba8-99c2b6d1f790
# ╠═01788aa6-486f-4bee-941c-c42a4e2a5f8e
# ╠═0a7010e7-c093-4c26-afd3-cf4c17af7049
# ╟─a0664415-813f-4e35-a519-ab44437720b9
# ╟─34ac14f9-752a-4095-866f-0f32befc6957
# ╠═46304beb-13f2-45d3-aef5-e9933472d42c
# ╠═3190dfb9-5bb3-41c9-81ae-1fcae71cd0a7
# ╠═b036d4ca-d8ea-43e1-ab23-bb99bbf91cec
# ╠═f0c29b93-c63c-47cf-bcbd-9b7b58cb7e99
# ╠═9b2379cd-7dff-4f5d-9c05-0e7ef4a707a0
# ╠═55e85624-6196-4501-b04a-f259aeb84752
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
