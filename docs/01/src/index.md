```@meta
EditURL = "https://github.com/TRAVIS_REPO_SLUG/blob/master/"
```

# Who am I ?

 - My name is *Pierre Navaro*
 - **Fortran 77 + PVM** : during my PhD 1998-2002 (Université du Havre)
 - **Fortran 90-2003 + OpenMP-MPI** : Engineer in Strasbourg (2003-2015) at IRMA
 - **Numpy + Cython, R + Rcpp** : engineer in Rennes (2015-now) at IRMAR
 - **Julia v1.0** since July 2018

---

# Why Julia?

- Born in 2009 and version 1.0 was released in August 2018.
- High-level languages like Python and R let one explore and experiment rapidly, but can run slow.
- Low-level languages like Fortran/C++ tend to take longer to develop, but run fast.
- This is sometimes called the "two language problem" and is something the Julia developers set out to eliminate.
- Julia's promise is to provide a "best of both worlds" experience for programmers who need to develop novel algorithms and bring them into production environments with minimal effort.

---

# Julia's Engineering and Design Tradoffs

- Type structures cannot bechanges after created (less dynamism but memory layout can be optimized for)
- All functions are JIT compiled via LLVM (interactive lags but massive runtime improvements)
- All functions specialize on types of arguments (more compilation but give generic programming structures)
- Julia is interactive (use it like Python and R, but makes it harder to get binaries)
- Julia has great methods for handling mutation (more optimization opportunities like C/Fortran, but more cognative burden)
- Julia's Base library and most packages are written in Julia, (you can understand the source, but learn a new package)
- Julia has expensive tooling for code generation and metaprogramming (consise and more optimizations, but some codes can be for experienced users)

To me, this gives me a language with a lot of depth which works well for computationally-expensive scientific applications.

[© ChrisRackaukas](https://www.youtube.com/watch?v=zJ3R6vOhibA&feature=em-uploademail)

---

# Type-Dispatch Programming

- Centered around implementing the generic template of the algorithm not around building representations of data.
- The data type choose how to efficiently implement the algorithm.
- With this feature share and reuse code is very easy

[JuliaCon 2019 | The Unreasonable Effectiveness of Multiple Dispatch | Stefan Karpinski](https://youtu.be/kc9HwsxE1OY)

---

```@example index
using DifferentialEquations, Plots

g = 9.79 # Gravitational constants
L = 1.00 # Length of the pendulum

#Initial Conditions
u₀ = [0, π / 60] # Initial speed and initial angle
tspan = (0.0, 6.3) # time domain

#Define the problem
function simplependulum(du, u, p, t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g/L)*θ
end

#Pass to solvers
prob = ODEProblem(simplependulum, u₀, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-6)
nothing # hide
```

---

Analytic solution

```@example index
u = u₀[2] .* cos.(sqrt(g / L) .* sol.t)

plot(sol.t, getindex.(sol.u, 2), label = "Numerical")
scatter!(sol.t, u, label = "Analytic")
savefig("pendulum1.svg")
```

![](pendulum1.svg)

---

[Numbers with Uncertainties](http://tutorials.juliadiffeq.org/html/type_handling/02-uncertainties.html)

```@example index
using Measurements

g = 9.79 ± 0.02; # Gravitational constants
L = 1.00 ± 0.01; # Length of the pendulum

#Initial Conditions
u₀ = [0 ± 0, π / 60 ± 0.01] # Initial speed and initial angle

#Define the problem
function simplependulum(du, u, p, t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g/L)*θ
end

#Pass to solvers
prob = ODEProblem(simplependulum, u₀, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-6);
nothing # hide
```

---

Analytic solution

```@example index
u = u₀[2] .* cos.(sqrt(g / L) .* sol.t)

plot(sol.t, getindex.(sol.u, 2), label = "Numerical")
plot!(sol.t, u, label = "Analytic")
savefig("pendulum2.svg");
```

![](pendulum2.svg)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

