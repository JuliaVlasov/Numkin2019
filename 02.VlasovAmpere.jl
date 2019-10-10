# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.2.0
#     language: julia
#     name: julia-1.2
# ---

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1D1V Vlasov–Ampere system
#
# $$
# \frac{\partial f}{\partial t} + \upsilon \frac{\partial f}{\partial x}
# - E(t,x) \frac{\partial f}{\partial \upsilon} = 0
# $$
#
# $$
# \frac{\partial E}{\partial t} = - J = \int f\upsilon \; d\upsilon
# $$

# + {"slideshow": {"slide_type": "slide"}}
using ProgressMeter, FFTW, Plots, LinearAlgebra
using BenchmarkTools, Statistics

# + {"slideshow": {"slide_type": "slide"}}
"""
    UniformMesh(start, stop, length)

1D uniform mesh data. This is a periodic domain so the end point is 
removed
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


# + {"slideshow": {"slide_type": "slide"}}
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

# + {"slideshow": {"slide_type": "slide"}}
"""
compute electric field from ρ
"""
function compute_e(mesh::UniformMesh, ρ)

   n = mesh.length
   k =  2π / (mesh.stop - mesh.start)
   modes = [1.0; k .* vcat(1:n÷2-1,-n÷2:-1)...]
   ρ̂ = fft(ρ)./modes
   vec(real(ifft(-1im .* ρ̂)))

end


# + {"slideshow": {"slide_type": "slide"}}
"""
    advection! = AmpereAdvection( mesh, kx)

    ∂f/∂t − v ∂f/∂x  = 0
    ∂E/∂t = −J = ∫ fv dv
    ∂f/∂t − E(x) ∂f/∂v  = 0

For every created struct, two methods are available
- Advection method along v
- Advection method along x and e computation

## Algorithm to compute electric field during advection along x

- For each ``j`` compute discrete Fourier transform in ``x``
 of ``f^n(x_i,\\upsilon_j)`` yielding ``f_k^n(\\upsilon_j)``, 
- For `` k \\neq 0 ``, compute 
```math
f^{n+1}_k(\\upsilon_j) = e^{−2i\\pi k \\upsilon
    \\Delta t/L} f_n^k(\\upsilon_j),
```
```math
\\rho_k^{n+1} = \\Delta \\upsilon \\sum_j􏰄 f^{n+1}_k(\\upsilon_j),
```
```math
E^{n+1}_k = \\rho^{n+1}_k L/(2i\\pi k \\epsilon_0),
```
- For ``k = 0`` do nothing: 
```math
f_{n+1}(\\upsilon_j) = f^n_k(\\upsilon_j), E^{n+1}_k = E^n_k.
```
- Perform inverse discrete Fourier transform of ``E^{n+1}_k`` and for each ``j`` of ``f^{n+1}_k (\\upsilon_j)``.
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


# + {"slideshow": {"slide_type": "slide"}}
?AmpereAdvection

# + {"slideshow": {"slide_type": "slide"}}

function (adv :: AmpereAdvection)( fᵗ  :: Array{ComplexF64,2}, 
                                   e   :: Vector{ComplexF64}, 
                                   dt  :: Float64 )
    fft!(fᵗ, 1)
    fᵗ .= fᵗ .* exp.(-1im * dt * adv.kx * transpose(e))
    ifft!(fᵗ, 1)

end


# + {"slideshow": {"slide_type": "slide"}}

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
        e[i] = -1im * ρ[i] ./ adv.kx[i]
    end
    e[1] = 0.0
    ifft!(f,1)
    ifft!(e)
    e .= real(e)
end


# + {"slideshow": {"slide_type": "slide"}}
"""
Landau damping initialisation function

[Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)

```math
f(x,v) = \\frac{1}{\\sqrt{2\\pi}}(1+ ϵ \\cdot cos(kₓ x)) e^{-v^2/2}
```

"""
function landau( ϵ, kx, x, v )
    
    (1.0.+ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.*v))
    
end

# + {"slideshow": {"slide_type": "slide"}}
?landau

# + {"slideshow": {"slide_type": "slide"}}
nx, nv = 256, 256
xmin, xmax =  0., 4*π
vmin, vmax = -6., 6.
tf = 60
nt = 600
meshx = UniformMesh(xmin, xmax, nx)
meshv = UniformMesh(vmin, vmax, nv)
            
# Initialize distribution function
x = meshx.points
v = meshv.points
ϵ, kx = 0.001, 0.5
    
# Allocate arrays for distribution function and its transposed
f = zeros(Complex{Float64},(nx,nv))
fᵀ= zeros(Complex{Float64},(nv,nx))
    
f .= landau( ϵ, kx, x, v)
    
transpose!(fᵀ,f)
    
ρ  = compute_rho(meshv, f)
e  = zeros(ComplexF64, nx)
e .= compute_e(meshx, ρ)
    
mod_e = Float64[]
    
dt = tf / nt
    
advection_x! = AmpereAdvection( meshx )
advection_v! = AmpereAdvection( meshv );

# +
        
@showprogress 1 for i in 1:nt
    
    advection_v!(fᵀ, e,  0.5dt)
    
    transpose!(f,fᵀ)
    
    advection_x!( f, e, v, dt)
    
    push!(mod_e, log(sqrt((sum(e.^2))*meshx.step)))
    
    transpose!(fᵀ,f)
    
    advection_v!(fᵀ, e,  0.5dt)
end

# + {"slideshow": {"slide_type": "slide"}}
t =  range(0,stop=tf,length=nt)
plot(t, -0.1533*t.-5.48)
plot!(t, mod_e , label=:ampere )
