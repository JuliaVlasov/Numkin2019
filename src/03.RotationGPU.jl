# ## GPU Computing

# https://github.com/JuliaGPU

using Plots, BenchmarkTools, FFTW, LinearAlgebra

# ### Advection equation for a rotation in two dimensional domain
#
# ```math
#  \\frac{d f}{dt} +  (y \\frac{d f}{dx} - x \\frac{d f}{dy}) = 0
# ```
#
# ``x \in [-π, π], y \in [-π, π] `` and  `` t \in [0, 200π] ``

#md # ---

struct Mesh
    
    nx   :: Int64
    ny   :: Int64
    x    :: Vector{Float64}
    y    :: Vector{Float64}
    kx   :: Vector{Float64}
    ky   :: Vector{Float64}
    
    function Mesh( xmin, xmax, nx, ymin, ymax, ny)
        ## periodic boundary condition, we remove the end point.
        x = range(xmin, stop=xmax, length=nx+1)[1:end-1]  
        y = range(ymin, stop=ymax, length=ny+1)[1:end-1]  
        kx  = 2π ./ (xmax-xmin) .* [0:nx÷2-1;nx÷2-nx:-1]
        ky  = 2π ./ (ymax-ymin) .* [0:ny÷2-1;ny÷2-ny:-1]
        new( nx, ny, x, y, kx, ky)
    end
end

#md # ---

function exact(time, mesh :: Mesh; shift=1.0)
   
    f = zeros(Float64,(mesh.nx, mesh.ny))
    for (i, x) in enumerate(mesh.x), (j, y) in enumerate(mesh.y)
        xn = cos(time)*x - sin(time)*y
        yn = sin(time)*x + cos(time)*y
        f[i,j] = exp(-(xn-shift)*(xn-shift)/0.1)*exp(-(yn-shift)*(yn-shift)/0.1)
    end

    f
end

#md # ---
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
       plot!(p[1]; clims=(0.,1.), aspect_ratio=:equal, colorbar=false, show=false)
       plot!(sqrt(2) .* cos.(-pi:0.1:pi+0.1), 
             sqrt(2) .* sin.(-pi:0.1:pi+0.1), label="", show=false)
       xlims!(-π,π)
       ylims!(-π,π)
       next!(bar) ## increment the progress bar
        
    end

    anim
    
end

#md # ---

anim = animation( 2π, 100)
#md gif(anim, "rotation2d.gif", fps = 30)
#md nothing # hide

#md # ![](rotation2d.gif)

#md # ---

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

#md # ---

mesh = Mesh( -π, π, 1024, -π, π, 1024)
nt, tf = 100, 20.
rotation_on_cpu(mesh, 1, 0.1)
@time norm( rotation_on_cpu(mesh, nt, tf) .- exact( tf, mesh))

#md # ---

using Pkg 

GPU_ENABLED = haskey(Pkg.installed(), "CUDAdrv")

if GPU_ENABLED

    using CUDAdrv, CuArrays, CuArrays.CUFFT
    
    println(CUDAdrv.name(CuDevice(0)))

end

#md # ---
    
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

#md # ---

if GPU_ENABLED

    nt, tf = 100, 20.
    rotation_on_gpu(mesh, 1, 0.1)
    @time norm( rotation_on_gpu(mesh, nt, tf) .- exact( tf, mesh))

end
