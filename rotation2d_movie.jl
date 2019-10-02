# -*- coding: utf-8 -*-
# # Who am I ?
#
# - My name is *Pierre Navaro*
# - Ph.D in Computational Aeroacoustics, 1998-2002 (Université du Havre) (Fortran 77+PVM)
# - Scientific Software Engineer in Strasbourg (2003-2015) (Fortran 90-2003 + OpenMP-MPI)
# - Moved to Rennes in 2015 (Numpy + Cython, R + Rcpp)
# - Julia user since July 2018 (Julia v1.0)
#     * Simulation numérique en physique des plasmas
#     * Assimilation de données météo

# # Advection equation for a rotation in two dimensional domain
#
# $$
# \frac{d f}{dt} +  (y \frac{d f}{dx} - x \frac{d f}{dy}) = 0
# $$
#
# $$ 
# x \in [-\pi, \pi],\qquad y \in [-\pi, \pi] \qquad \mbox{ and } \qquad t \in [0, 200\pi] 
# $$
#
#
# https://github.com/JuliaVlasov/FourierAdvections.jl in `notebooks` directory
#
# You can open these julia files as notebooks with [jupytext](https://jupytext.readthedocs.io)

using FFTW, LinearAlgebra, Plots, ProgressMeter
using BenchmarkTools

"""
    Mesh( xmin, xmax, nx, ymin, ymax, ny)

mesh information
"""
struct Mesh
    
    nx   :: Int64
    ny   :: Int64
    x    :: Vector{Float64}
    y    :: Vector{Float64}
    kx   :: Vector{Float64}
    ky   :: Vector{Float64}
    
    function Mesh( xmin, xmax, nx, ymin, ymax, ny)
        
        x = range(xmin, stop=xmax, length=nx+1)[1:end-1]  ## we remove the end point
        y = range(ymin, stop=ymax, length=ny+1)[1:end-1]  ## periodic boundary condition
        kx = 2π/(xmax-xmin)*[0:nx÷2-1;nx÷2-nx:-1]
        ky = 2π/(ymax-ymin)*[0:ny÷2-1;ny÷2-ny:-1]

        new( nx, ny, x, y, kx, ky)
    end
end

# ### Function to compute exact solution

function exact(time, mesh :: Mesh; shift=1)
   
    f = zeros(Float64,(mesh.nx, mesh.ny))
    for (i, x) in enumerate(mesh.x), (j, y) in enumerate(mesh.y)  # two loops
        xn = cos(time)*x - sin(time)*y
        yn = sin(time)*x + cos(time)*y
        f[i,j] = exp(-(xn-shift)*(xn-shift)/0.1)*exp(-(yn-shift)*(yn-shift)/0.1)
    end

    f
end

mesh = Mesh(-π, π, 128, -π, π, 128)
f = exact(0.0, mesh)
contour(mesh.x, mesh.y, f; aspect_ratio=:equal, clims=(0.,1.))

# ## Create the gif to show what we are computing

# +
function animation( tf, nt)
    
    nx, ny = 64, 64
    xmin, xmax, nx = -π, π, nx
    ymin, ymax, ny = -π, π, ny
    mesh = Mesh(xmin, xmax, nx, ymin, ymax, ny)
    f  = zeros(Float64,(nx,ny))
    dt = tf / nt
    bar = Progress(nt,1) ## progress bar
    t = 0
    anim = @animate for n=1:nt
       
       f .= exact(t, mesh)
       t += dt
       p = contour(mesh.x, mesh.y, f)
       plot!(p[1]; clims=(0.,1.), aspect_ratio=:equal, colorbar=false)
       plot!(sqrt(2) .* cos.(-pi:0.1:pi+0.1), sqrt(2) .* sin.(-pi:0.1:pi+0.1))
       next!(bar) ## increment the progress bar
        
    end every 1

    anim
    
end

anim = animation( 2π, 100)
gif(anim, "rotation2d.gif", fps = 30)
