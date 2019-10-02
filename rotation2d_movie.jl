using FFTW, LinearAlgebra, ProgressMeter, BenchmarkTools

using Plots
default(show=false)


# ### Function to compute exact solution

function exact!(f, time, x, y; shift=1)
   
    nx = length(x)
    ny = length(y)
    for (i, x) in enumerate(x), (j, y) in enumerate(y)
        xn = cos(time)*x - sin(time)*y
        yn = sin(time)*x + cos(time)*y
        f[i,j] = (exp(-(xn-shift)*(xn-shift)/0.1)
                 *exp(-(yn-shift)*(yn-shift)/0.1))
    end

end


function animation( tf, nt)
    
    nx, ny = 64, 64
    xmin, xmax, nx = -π, π, nx
    ymin, ymax, ny = -π, π, ny
    x = range(xmin, stop=xmax, length=nx+1)[1:end-1] 
    y = range(ymin, stop=ymax, length=ny+1)[1:end-1]
    f  = zeros(Float64,(nx,ny))
    dt = tf / nt
    bar = Progress(nt,1) ## progress bar
    t = 0
    anim = @animate for n=1:nt
       
       exact!(f, t, x, y)
       t += dt
       p = contour(x, y, f, axis=[], framestyle=:none)
       plot!(p[1]; clims=(0.,1.), aspect_ratio=:equal, colorbar=false)
       plot!(sqrt(2) .* cos.(-pi:0.1:pi+0.1), 
             sqrt(2) .* sin.(-pi:0.1:pi+0.1), label="")
       xlims!(-π,π)
       ylims!(-π,π)
       next!(bar) ## increment the progress bar
        
    end

    anim
    
end

anim = animation( 2π, 100)
gif(anim, "rotation2d.gif", fps = 30)
