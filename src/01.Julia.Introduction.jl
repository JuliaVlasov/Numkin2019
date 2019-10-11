#md # # Who am I ?
#md #
#md #  - My name is *Pierre Navaro*
#md #  - **Fortran 77 + PVM** : during my PhD 1998-2002 (Université du Havre) 
#md #  - **Fortran 90-2003 + OpenMP-MPI** : Engineer in Strasbourg (2003-2015) at IRMA
#md #  - **Numpy + Cython, R + Rcpp** : engineer in Rennes (2015-now) at IRMAR
#md #  - **Julia v1.0** since July 2018 

#md # ---

#md #
#md # # Why Julia?
#md #
#md # - Born in 2009 and version 1.0 was released in August 2018.
#md # - High-level languages like Python and R let one explore and experiment rapidly, but can run slow.
#md # - Low-level languages like Fortran/C++ tend to take longer to develop, but run fast.
#md # - This is sometimes called the "two language problem" and is something the Julia developers set out to eliminate. 
#md # - Julia's promise is to provide a "best of both worlds" experience for programmers who need to develop novel algorithms and bring them into production environments with minimal effort.
#md #

#md # ---

#md # # Julia's Engineering and Design Tradoffs
#md #
#md # - Type structures cannot bechanges after created (less dynamism but memory layout can be optimized for)
#md # - All functions are JIT compiled via LLVM (interactive lags but massive runtime improvements)
#md # - All functions specialize on types of arguments (more compilation but give generic programming structures)
#md # - Julia is interactive (use it like Python and R, but makes it harder to get binaries)
#md # - Julia has great methods for handling mutation (more optimization opportunities like C/Fortran, but more cognative burden)
#md # - Julia's Base library and most packages are written in Julia, (you can understand the source, but learn a new package)
#md # - Julia has expensive tooling for code generation and metaprogramming (consise and more optimizations, but some codes can be for experienced users)
#md #
#md # To me, this gives me a language with a lot of depth which works well for computationally-expensive scientific applications.
#md #
#md # [© ChrisRackaukas](https://www.youtube.com/watch?v=zJ3R6vOhibA&feature=em-uploademail) 

#md # ---

#md # # Type-Dispatch Programming
#md #
#md # - Centered around implementing the generic template of the algorithm not around building representations of data.
#md # - The data type choose how to efficiently implement the algorithm.
#md # - With this feature share and reuse code is very easy
#md #
#md # [JuliaCon 2019 | The Unreasonable Effectiveness of Multiple Dispatch | Stefan Karpinski](https://youtu.be/kc9HwsxE1OY)

#md # ---

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

#md # ---

# Analytic solution
u = u₀[2] .* cos.(sqrt(g / L) .* sol.t)

plot(sol.t, getindex.(sol.u, 2), label = "Numerical")
scatter!(sol.t, u, label = "Analytic")
#md savefig("pendulum1.svg")

#md # ![](pendulum1.svg)

#md # ---

# [Numbers with Uncertainties](http://tutorials.juliadiffeq.org/html/type_handling/02-uncertainties.html)

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

#md # ---

# Analytic solution
u = u₀[2] .* cos.(sqrt(g / L) .* sol.t)

plot(sol.t, getindex.(sol.u, 2), label = "Numerical")
plot!(sol.t, u, label = "Analytic")
#md savefig("pendulum2.svg"); 


#md # ![](pendulum2.svg)


#md # ---
