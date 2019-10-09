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

using Plots, BenchmarkTools, FFTW, LinearAlgebra

# ### Advection equation for a rotation in two dimensional domain
#
#  $$
#  \frac{d f}{dt} +  (y \frac{d f}{dx} - x \frac{d f}{dy}) = 0
#  $$
#
#  $$
#  x \in [-\pi, \pi],\qquad y \in [-\pi, \pi] \qquad \mbox{ and } \qquad t \in [0, 200\pi]
#  $$

# ![](rotation2d.gif)

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
        kx  = 2π ./ (xmax-xmin) .* [0:nx÷2-1;nx÷2-nx:-1]
        ky  = 2π ./ (ymax-ymin) .* [0:ny÷2-1;ny÷2-ny:-1]
        new( nx, ny, x, y, kx, ky)
    end
end

function exact(time, mesh :: Mesh; shift=1.0)
   
    f = zeros(Float64,(mesh.nx, mesh.ny))
    for (i, x) in enumerate(mesh.x), (j, y) in enumerate(mesh.y)
        xn = cos(time)*x - sin(time)*y
        yn = sin(time)*x + cos(time)*y
        f[i,j] = exp(-(xn-shift)*(xn-shift)/0.1)*exp(-(yn-shift)*(yn-shift)/0.1)
    end

    f
end

# +
function rotation_on_cpu( mesh :: Mesh, nt :: Int64, tf :: Float64) 
    
    dt = tf / nt
    
    f   = zeros(ComplexF64,(mesh.nx,mesh.ny))
    f  .= exact( 0.0, mesh )
    
    exky = exp.( 1im*tan(dt/2) .* mesh.x  .* mesh.ky')
    ekxy = exp.(-1im*sin(dt)   .* mesh.y' .* mesh.kx )
    
    for n = 1:nt
        
        fft!(f, 2)
        f .= exky .* f  # 
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
# -

mesh = Mesh( -π, π, 1024, -π, π, 1024)
nt, tf = 100, 20.
rotation_on_cpu(mesh, 1, 0.1)
@time norm( rotation_on_cpu(mesh, nt, tf) .- exact( tf, mesh))

using CUDAdrv
CUDAdrv.name(CuDevice(0))

using CuArrays
using CuArrays.CUFFT

# +
function rotation_on_gpu( mesh :: Mesh, nt :: Int64, tf :: Float64)
    
    dt = tf / nt
    
    f   = zeros(ComplexF64,(mesh.nx, mesh.ny))
    f  .= exact( 0.0, mesh)
    
    d_f    = CuArray(f) # allocate f on GPU
    
    # Create fft plans on GPU
    p_x    = plan_fft!(d_f,  [1])
    pinv_x = plan_ifft!(d_f, [1])
    p_y    = plan_fft!(d_f,  [2])
    pinv_y = plan_ifft!(d_f, [2])  
    
    # Allocate on GPU
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
# -

nt, tf = 100, 20.
rotation_on_gpu(mesh, 1, 0.1)
@time norm( rotation_on_gpu(mesh, nt, tf) .- exact( tf, mesh))
