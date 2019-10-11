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
## 1D1V Vlasov–Ampere system

```math
\\frac{\\partial f}{\\partial t} + \\upsilon \\frac{\\partial f}{\\partial x} - E(t,x) \\frac{\\partial f}{\\partial \\upsilon} = 0
```

```math
\\frac{\\partial E}{\\partial t} = - J = \\displaystyle \\int  f\\upsilon \\; d\\upsilon
```

---

```@example index
using ProgressMeter, FFTW, Plots, LinearAlgebra
using BenchmarkTools, Statistics

"""
    UniformMesh(start, stop, length)

1D uniform mesh data for periodic domain (end point is removed)
"""
struct UniformMesh

   start    :: Float64
   stop     :: Float64
   length   :: Int64
   step     :: Float64
   points   :: Vector{Float64}

   function UniformMesh(start, stop, length)

       points = range(start, stop=stop, length=length+1)[1:end-1]
       step = (stop - start) / length

       new( start, stop, length, step, points)

   end

end
```

---

## Compute charge density ρ(x)

```@example index
"""
    compute_rho( mesh, f)

Compute charge density ρ(x,t) = ∫ f(x,v,t) dv

returns ρ - ρ̄

"""
function compute_rho(meshv::UniformMesh, f)

   dv  = meshv.step
   ρ = dv .* vec(sum(real(f), dims=2))
   ρ .- mean(ρ)

end
```

---

## Compute electric field from ρ(x)

```@example index
"""
    compute_e(mesh, ρ)

compute electric field from ρ
"""
function compute_e(mesh::UniformMesh, ρ)

   n = mesh.length
   k =  2π / (mesh.stop - mesh.start)
   modes = [1.0; k .* vcat(1:n÷2-1,-n÷2:-1)...]
   ρ̂ = fft(ρ)./modes
   vec(real(ifft(-1im .* ρ̂)))

end
```

---

## Callable struct `Advection`

```@example index
"""
    advection! = AmpereAdvection( mesh, kx)

For every created struct, two methods are available
- Advection method along v
- Advection method along x and e computation

"""
struct AmpereAdvection

    mesh :: UniformMesh
    kx   :: Vector{Float64}

    function AmpereAdvection( mesh )

        nx  = mesh.length
        dx  = mesh.step
        Lx  = mesh.stop - mesh.start
        kx  = zeros(Float64, nx)
        kx .= 2π/Lx .* [0:nx÷2-1;-nx÷2:-1]
        new( mesh, kx)

    end

end
```

---

```@example index
function (adv :: AmpereAdvection)( fᵗ  :: Array{ComplexF64,2},
                                   e   :: Vector{ComplexF64},
                                   dt  :: Float64 )
    fft!(fᵗ, 1)
    fᵗ .= fᵗ .* exp.(-1im * dt * adv.kx * transpose(e))
    ifft!(fᵗ, 1)

end
```

4UkInfedi7B7lKAwjXdpLEd52ajQ8FCDMXGjPBurjZQXb6xR6V

```@example index
function (adv :: AmpereAdvection)( f   :: Array{ComplexF64,2},
                                   e   :: Vector{ComplexF64},
                                   v   :: Vector{Float64},
                                   dt  :: Float64 )

    ev = exp.(-1im*dt * adv.kx * transpose(v))

    fft!(f,1)
    f .= f .* ev
    dv = v[2]-v[1]
    ρ = dv * vec(sum(f,dims=2))
    for i in 2:length(e)
        e[i] = -1im * ρ[i] / adv.kx[i]
    end
    e[1] = 0.0
    ifft!(f, 1)
    ifft!(e)
    e .= real(e)
end
```

----

## Initial distribution function

```math
f(x,v) = \\frac{1}{\\sqrt{2\\pi}}(1+ ϵ \\cdot cos(k_x x)) e^{-v^2/2}
```

```@example index
"""
    landau( ϵ, kx, x, v )

Landau damping initialisation function

[Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)


"""
function landau( ϵ, kx, x, v )

    (1.0 .+ ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.*v))

end
```

---

```@example index
nx, nv = 256, 256
xmin, xmax =  0., 4*π
vmin, vmax = -6., 6.
tf = 60
nt = 600
meshx = UniformMesh(xmin, xmax, nx)
meshv = UniformMesh(vmin, vmax, nv);
```

4UkInfedi7B7lKAwjXdpLEd52ajQ8FCDMXGjPBurjZQXb6xR6V

Initialize distribution function

```@example index
x = meshx.points
v = meshv.points
ϵ, kx = 0.001, 0.5
```

Allocate arrays for distribution function and its transposed

```@example index
f = zeros(Complex{Float64},(nx,nv))
fᵀ= zeros(Complex{Float64},(nv,nx))

f .= landau( ϵ, kx, x, v);
```

---

Plot the distribution

```@example index
surface( x, v, real(f))
savefig("landau.svg")
```

![](landau.svg)

---

```@example index
transpose!(fᵀ,f)

ρ  = compute_rho(meshv, f)
e  = zeros(ComplexF64, nx)
e .= compute_e(meshx, ρ)

mod_e = Float64[]

dt = tf / nt

advection_x! = AmpereAdvection( meshx )
advection_v! = AmpereAdvection( meshv );
```

---

```@example index
@showprogress 1 for i in 1:nt

    advection_v!(fᵀ, e, 0.5dt)

    transpose!(f, fᵀ)

    advection_x!( f, e, v, dt)

    push!(mod_e, log(sqrt((sum(e.^2))*meshx.step))) # diagnostic

    transpose!(fᵀ, f)

    advection_v!(fᵀ, e, 0.5dt)

end
```

---

```@example index
t = range(0, stop=tf, length=nt)
plot(t, -0.1533*t.-5.48)
plot!(t, mod_e , label=:ampere )
savefig("mod_e.svg")
```

![](mod_e.svg)

---
## GPU Computing

https://github.com/JuliaGPU

```@example index
using Plots, BenchmarkTools, FFTW, LinearAlgebra
```

### Advection equation for a rotation in two dimensional domain

```math
 \\frac{d f}{dt} +  (y \\frac{d f}{dx} - x \\frac{d f}{dy}) = 0
```

``x \in [-π, π], y \in [-π, π] `` and  `` t \in [0, 200π] ``

---

```@example index
struct Mesh

    nx   :: Int64
    ny   :: Int64
    x    :: Vector{Float64}
    y    :: Vector{Float64}
    kx   :: Vector{Float64}
    ky   :: Vector{Float64}

    function Mesh( xmin, xmax, nx, ymin, ymax, ny)
        # periodic boundary condition, we remove the end point.
        x = range(xmin, stop=xmax, length=nx+1)[1:end-1]
        y = range(ymin, stop=ymax, length=ny+1)[1:end-1]
        kx  = 2π ./ (xmax-xmin) .* [0:nx÷2-1;nx÷2-nx:-1]
        ky  = 2π ./ (ymax-ymin) .* [0:ny÷2-1;ny÷2-ny:-1]
        new( nx, ny, x, y, kx, ky)
    end
end
```

---

```@example index
function exact(time, mesh :: Mesh; shift=1.0)

    f = zeros(Float64,(mesh.nx, mesh.ny))
    for (i, x) in enumerate(mesh.x), (j, y) in enumerate(mesh.y)
        xn = cos(time)*x - sin(time)*y
        yn = sin(time)*x + cos(time)*y
        f[i,j] = exp(-(xn-shift)*(xn-shift)/0.1)*exp(-(yn-shift)*(yn-shift)/0.1)
    end

    f
end
```

---

```@example index
using ProgressMeter

function animation( tf, nt)

    mesh = Mesh( -π, π, 64, -π, π, 64)
    dt = tf / nt
    bar = Progress(nt,1) ## progress bar
    t = 0
    anim = @animate for n=1:nt

       f = exact(t, mesh)
       t += dt
       p = contour(mesh.x, mesh.y, f, axis=[], framestyle=:none)
       plot!(p[1]; clims=(0.,1.), aspect_ratio=:equal, colorbar=false)
       plot!(sqrt(2) .* cos.(-pi:0.1:pi+0.1),
             sqrt(2) .* sin.(-pi:0.1:pi+0.1), label="")
       xlims!(-π,π)
       ylims!(-π,π)
       next!(bar) ## increment the progress bar

    end

    anim

end
```

---

```@example index
anim = animation( 2π, 100)
gif(anim, "rotation2d.gif", fps = 30)
nothing # hide
```

![](rotation2d.gif)

---

```@example index
function rotation_on_cpu( mesh :: Mesh, nt :: Int64, tf :: Float64)

    dt = tf / nt

    f   = zeros(ComplexF64,(mesh.nx,mesh.ny))
    f  .= exact( 0.0, mesh )

    exky = exp.( 1im*tan(dt/2) .* mesh.x  .* mesh.ky')
    ekxy = exp.(-1im*sin(dt)   .* mesh.y' .* mesh.kx )

    for n = 1:nt

        fft!(f, 2)
        f .= exky .* f
        ifft!(f,2)

        fft!(f, 1)
        f .= ekxy .* f
        ifft!(f, 1)

        fft!(f, 2)
        f .= exky .* f
        ifft!(f, 2)

    end

    real(f)

end
```

---

```@example index
mesh = Mesh( -π, π, 1024, -π, π, 1024)
nt, tf = 100, 20.
rotation_on_cpu(mesh, 1, 0.1)
@time norm( rotation_on_cpu(mesh, nt, tf) .- exact( tf, mesh))
```

---

```@example index
using Pkg

GPU_ENABLED = haskey(Pkg.installed(), "CUDAdrv")

if GPU_ENABLED

    using CUDAdrv, CuArrays, CuArrays.CUFFT

    println(CUDAdrv.name(CuDevice(0)))

end
```

---

```@example index
if GPU_ENABLED

    function rotation_on_gpu( mesh :: Mesh, nt :: Int64, tf :: Float64)

        dt = tf / nt

        f   = zeros(ComplexF64,(mesh.nx, mesh.ny))
        f  .= exact( 0.0, mesh)

        d_f    = CuArray(f) # allocate f on GPU

        p_x    = plan_fft!(d_f,  [1]) # Create fft plans on GPU
        pinv_x = plan_ifft!(d_f, [1])
        p_y    = plan_fft!(d_f,  [2])
        pinv_y = plan_ifft!(d_f, [2])

        d_exky = CuArray(exp.( 1im*tan(dt/2) .* mesh.x  .* mesh.ky'))
        d_ekxy = CuArray(exp.(-1im*sin(dt)   .* mesh.y' .* mesh.kx ))

        for n = 1:nt

            p_y * d_f
            d_f .*= d_exky
            pinv_y * d_f

            p_x * d_f
            d_f .*= d_ekxy
            pinv_x * d_f

            p_y * d_f
            d_f .*= d_exky
            pinv_y * d_f

        end

        f .= collect(d_f) # Transfer f from GPU to CPU

        real(f)

    end

end
```

---

```@example index
if GPU_ENABLED

    nt, tf = 100, 20.
    rotation_on_gpu(mesh, 1, 0.1)
    @time norm( rotation_on_gpu(mesh, nt, tf) .- exact( tf, mesh))

end
```

# Metaprogramming

## The ParticleGroup example

```@example index
import Sobol
using Plots, LinearAlgebra
```

    ParticleGroup{D,V}(n_particles, charge, mass)

- `D` : number of dimension in physical space
- `V` : number of dimension in phase space
- `n` : number of particles

```@example index
mutable struct ParticleGroup{D,V}

    n_particles       :: Int64
    data              :: Array{Float64, 2}

    function ParticleGroup{D,V}(n) where {D, V}

        data = zeros( Float64, (D+V, n))
        new( n, data)
    end
end
```

---

Set position of ith particle of p to x

```@example index
@generated function set_x!( p :: ParticleGroup{D,V}, i, x :: Float64 ) where {D, V}

    :(p.data[1, i] = x)

end
```

4UkInfedi7B7lKAwjXdpLEd52ajQ8FCDMXGjPBurjZQXb6xR6V

Set position of ith particle of p to x when x is a vector

```@example index
@generated function set_x!( p :: ParticleGroup{D,V}, i, x :: Vector{Float64} ) where {D, V}

     :(for j in 1:$D p.data[j, i] = x[j] end)

end
```

---

Set velocity of ith particle of p to v

```@example index
@generated function set_v!( p :: ParticleGroup{D,V}, i, v :: Float64 ) where {D, V}

    :(p.data[$D+1, i] = v)

end
```

4UkInfedi7B7lKAwjXdpLEd52ajQ8FCDMXGjPBurjZQXb6xR6V

Set velocity of ith particle of p to v

```@example index
@generated function set_v!( p :: ParticleGroup{D,V}, i, v :: Vector{Float64} ) where {D, V}

    :(for j in 1:$V p.data[$D+j, i] = v[j] end)

end
```

---

Get position of ith particle of p

```@example index
@generated function get_x( p :: ParticleGroup{D,V}, i ) where {D, V}

     :(p.data[1:$D, i])

end
```

Get velocity of ith particle of p

```@example index
@generated function get_v( p :: ParticleGroup{D,V}, i ) where {D, V}

     :(p.data[$D+1:$D+$V, i])
end
```

---

Sampling from a probability distribution to initialize Landau damping

```@example index
function landau_sampling!( pg :: ParticleGroup{1,2}, alpha, kx )

    function newton(r)
        x0, x1 = 0.0, 1.0
        r *= 2π / kx
        while (abs(x1-x0) > 1e-12)
            p = x0 + alpha * sin( kx * x0) / kx
            f = 1 + alpha * cos( kx * x0)
            x0, x1 = x1, x0 - (p - r) / f
        end
        x1
    end

    s = Sobol.SobolSeq(2)
    nbpart = pg.n_particles

    for i=1:nbpart

        v = sqrt(-2 * log( (i-0.5)/nbpart))
        r1, r2 = Sobol.next!(s)
        θ = r1 * 2π
        set_x!(pg,  i, newton(r2))
        set_v!(pg,  i, [v * cos(θ), v * sin(θ)])
    end

end

#+
```

---

```@example index
n_particles = 10000
pg = ParticleGroup{1,2}( n_particles)
alpha, kx = 0.1, 0.5
landau_sampling!(pg, alpha, kx)
```

4UkInfedi7B7lKAwjXdpLEd52ajQ8FCDMXGjPBurjZQXb6xR6V

```@example index
#+
xp = vcat([get_x(pg, i) for i in 1:pg.n_particles]...)
vp = vcat([get_v(pg, i) for i in 1:pg.n_particles]'...)
```

---

```@example index
#+
pp = plot(layout=(3,1))
histogram!(pp[1,1], xp, normalize=true, bins = 100, lab=:x)
plot!(pp[1,1], x -> (1+alpha*cos(kx*x))/(2π/kx), 0., 2π/kx, lab="")
histogram!(pp[2,1], vp[:,1], normalize=true, bins = 100, lab=:vx)
plot!(pp[2,1], v -> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
histogram!(pp[3,1], vp[:,2], normalize=true, bins = 100, lab=:vy)
plot!(pp[3,1], v -> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
savefig("particles.svg"); nothing # hide
```

![](particles.svg)

---

```@example index
#+
histogram2d(vp[:,1], vp[:,2], normalize=true, bins=100)
savefig("hist2d.svg")
```

![](hist2d.svg)

---
# GEMPIC

### Geometric ElectroMagnetic Particle-In-Cell Methods

https://arxiv.org/abs/1609.03053

Michael Kraus, Katharina Kormann, Philip J. Morrison, Eric Sonnendrücker

Framework for Finite Element Particle-in-Cell methods based on
the discretization of the underlying Hamiltonian structure of the
Vlasov-Maxwell system.

Install the GEMPIC package

```@example index
using Pkg

Pkg.add(PackageSpec(url="https://github.com/juliavlasov/GEMPIC.jl"))

using ProgressMeter, Plots, GEMPIC
```

---

# Strong Landau Damping

The physical parameters

```@example index
kx, α = 0.5, 0.5
xmin, xmax = 0, 2π/kx
domain = [xmin, xmax, xmax - xmin]
```

4UkInfedi7B7lKAwjXdpLEd52ajQ8FCDMXGjPBurjZQXb6xR6V

The numerical parameters

```@example index
∆t = 0.05
nx = 32
n_particles = 100000
mesh = GEMPIC.Mesh( xmin, xmax, nx)
spline_degree = 3
```

---

Initialize particles

```@example index
mass, charge = 1.0, 1.0
particle_group = GEMPIC.ParticleGroup{1,2}( n_particles, mass, charge, 1)
sampler = LandauDamping( α, kx)

sample!(sampler, particle_group)
```

4UkInfedi7B7lKAwjXdpLEd52ajQ8FCDMXGjPBurjZQXb6xR6V

Particle-mesh coupling operators

```@example index; continued = true
kernel_smoother1 = ParticleMeshCoupling( domain, [nx], n_particles,
                                         spline_degree-1, :galerkin)

kernel_smoother0 = ParticleMeshCoupling( domain, [nx], n_particles,
                                         spline_degree, :galerkin)
```

Allocate electrostatic fields and Maxwell solver

```@example index
rho = zeros(Float64, nx)
efield_poisson = zeros(Float64, nx)

maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)
```

---

### Charge density

```@example index
xg = LinRange(xmin, xmax, nx)
sval = eval_uniform_periodic_spline_curve(spline_degree-1, rho)
plot( xg, sval, label="ρ")
savefig("rho.svg")
```

![](rho.svg)

---

### Electric field

```@example index
solve_poisson!( efield_poisson, particle_group,
                kernel_smoother0, maxwell_solver, rho)

sval = eval_uniform_periodic_spline_curve(spline_degree-1, efield_poisson)
plot( xg, sval )
savefig("ex.svg")
```

![](ex.svg)

---

Initialize the arrays for the spline coefficients of the fields

```@example index; continued = true
efield_dofs = [copy(efield_poisson), zeros(Float64, nx)]
bfield_dofs = zeros(Float64, nx)

propagator = HamiltonianSplitting( maxwell_solver,
                                   kernel_smoother0,
                                   kernel_smoother1,
                                   particle_group,
                                   efield_dofs,
                                   bfield_dofs,
                                   domain);

efield_dofs_n = propagator.e_dofs

thdiag = TimeHistoryDiagnostics( particle_group, maxwell_solver,
                        kernel_smoother0, kernel_smoother1 );
```

---

## Loop over time

```@example index
steps, Δt = 100, 0.05

@showprogress 1 for j = 1:steps # loop over time

    # Strang splitting
    strang_splitting!(propagator, Δt, 1)

    # Diagnostics
    solve_poisson!( efield_poisson, particle_group,
                    kernel_smoother0, maxwell_solver, rho)

    write_step!(thdiag, j * Δt, spline_degree,
                    efield_dofs,  bfield_dofs,
                    efield_dofs_n, efield_poisson)

end
```

---

## Diagnostics stored in a dataframe

```@example index
using DataFrames
first(thdiag.data, 5)
```

---

```@example index
import Gadfly: Geom, Scale

Gadfly.plot(thdiag.data, x=:Time, y=:PotentialEnergyE1, Geom.line, Scale.y_log10)
savefig("thdiag.svg")
nothing # hide
```

![](thdiag.svg)

---
### Optimizing Julia code is often done at the expense of transparency

```@example index
using Random, LinearAlgebra, BenchmarkTools

function test(A, B, C)
    C = C - A * B
    return C
end

A = rand(1024, 256); B = rand(256, 1024); C = rand(1024, 1024)
@btime test(A, B, C); #C, A and B are matrices.
nothing # hide
```

4UkInfedi7B7lKAwjXdpLEd52ajQ8FCDMXGjPBurjZQXb6xR6V

```@example index
function test_opt(A, B, C)
    BLAS.gemm!('N','N', -1., A, B, 1., C)
    return C
end
@btime test_opt(A, B, C) # avoids taking two unnecessary copies of the matrix C.
nothing # hide
```

4UkInfedi7B7lKAwjXdpLEd52ajQ8FCDMXGjPBurjZQXb6xR6V

```@example index
C = rand(1024, 1024)
all(test(A, B, C) .== test_opt(A, B, C))
```

---

### Derivative computation with FFT

```@example index
using FFTW

xmin, xmax, nx = 0, 4π, 1024
ymin, ymax, ny = 0, 4π, 1024

x = range(xmin, stop=xmax, length=nx+1)[1:end-1]
y = range(ymin, stop=ymax, length=ny+1)[1:end-1]

ky  = 2π ./ (ymax-ymin) .* [0:ny÷2-1;ny÷2-ny:-1]
exky = exp.( 1im .* ky' .* x)

function df_dy( f, exky )
    ifft(exky .* fft(f, 2), 2)
end

f = sin.(x) .* cos.(y') # f is a 2d array created by broadcasting

@btime df_dy(f, exky)
nothing # hide
```

---

### Memory alignement, and inplace computation.

```@example index
f  = zeros(ComplexF64, (nx,ny))
fᵗ = zeros(ComplexF64, reverse(size(f)))
f̂ᵗ = zeros(ComplexF64, reverse(size(f)))

f .= sin.(x) .* cos.(y')

fft_plan = plan_fft(fᵗ, 1, flags=FFTW.PATIENT)

function df_dy!( f, fᵗ, f̂ᵗ, exky )
    transpose!(fᵗ,f)
    mul!(f̂ᵗ,  fft_plan, fᵗ)
    f̂ᵗ .= f̂ᵗ .* exky
    ldiv!(fᵗ, fft_plan, f̂ᵗ)
    transpose!(f, fᵗ)
end

@btime df_dy!(f, fᵗ, f̂ᵗ, exky )
nothing # hide
```

---

# Why use Julia language!

- **You develop in the same language in which you optimize.**
- Packaging system is very efficient (3173 registered packages)
- PyPi (198,360 projects) R (14993 packages)
- It is very easy to create a package (easier than R and Python)
- It is very easy to use GPU device.
- Nice interface for Linear Algebra and Differential Equations
- Easy access to BLAS and LAPACK
- Julia talks to all major Languages - mostly without overhead!

---

# What's bad

- It is still hard to build shared library or executable from Julia code.
- Compilation times can be unbearable.
- Plotting takes time (20 seconds for the first plot)
- OpenMP is better than the Julia multithreading library but it is progressing.
- There is a MPI and PETSc package but they are not very active.
- For parallelization, The Julia community seems to prefer the ditributed processing approach.
- Does not work well with vectorized code, you need to do a lot of inplace computation to avoid memory allocations and use explicit views to avoid copy.
- Julia website proclaims that it is faster than Fortran but this is not true.

[What's Bad About Julia by Jeff Bezanson](https://www.youtube.com/watch?v=TPuJsgyu87U)

---

## So when should i use Julia?

- Now! If you need performance and plan to write your own libraries.
- In ~1-2 Years if you want a smooth deploy.
- In ~3-5 Years if you want a 100% smooth experience.

## Julia Munich Meetup
Every two months, poll for the next meeting that
will take place at Garching Campus : https://doodle.com/poll/z3ft2dytnaebyhh7.

## Python-Julia benchmarks by Thierry Dumont

https://github.com/Thierry-Dumont/BenchmarksPythonJuliaAndCo/wiki

---

# Julia is a language made for Science.

 [Some State of the Art Packages](http://www.stochasticlifestyle.com/some-state-of-the-art-packages-in-julia-v1-0)

 * JuliaDiffEq – Differential equation solving and analysis.
 * JuliaDiff – Differentiation tools.
 * JuliaGeometry – Computational Geometry.
 * JuliaGraphs – Graph Theory and Implementation.
 * JuliaIntervals - Rigorous numerics with interval arithmetic & applications.
 * JuliaMath – Mathematics made easy in Julia.
 * JuliaOpt – Optimization.
 * JuliaPolyhedra – Polyhedral computation.
 * JuliaSparse – Sparse matrix solvers.
 * JuliaStats – Statistics and Machine Learning.
 * JuliaPlots - powerful convenience for visualization.
 * JuliaGPU - GPU Computing for Julia.
 * FluxML - The Elegant Machine Learning Stack.

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

