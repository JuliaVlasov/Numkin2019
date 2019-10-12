ENV["GKSwstype"]="100" #src

# ## 1D1V Vlasov–Ampere system
#
# ```math
# \\frac{\\partial f}{\\partial t} + \\upsilon \\frac{\\partial f}{\\partial x} - E(t,x) \\frac{\\partial f}{\\partial \\upsilon} = 0
# ```
#
# ```math
# \\frac{\\partial E}{\\partial t} = - J = \\displaystyle \\int  f\\upsilon \\; d\\upsilon
# ```

#md # ---

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

#md # ---
#
# ## Compute charge density ρ(x)
#

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

#md # ---
#
# ## Compute electric field from ρ(x)
#

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

#
#md # ---
#
# ## Callable struct `Advection`
#
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


#md # ---

function (adv :: AmpereAdvection)( fᵗ  :: Array{ComplexF64,2}, 
                                   e   :: Vector{ComplexF64}, 
                                   dt  :: Float64 )
    fft!(fᵗ, 1)
    fᵗ .= fᵗ .* exp.(-1im * dt * adv.kx * transpose(e))
    ifft!(fᵗ, 1)

end

#md # --

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

#
#
#md # ----
#
# ## Initial distribution function
#
# ```math
# f(x,v) = \\frac{1}{\\sqrt{2\\pi}}(1+ ϵ \\cdot cos(k_x x)) e^{-v^2/2}
# ```
#
"""
    landau( ϵ, kx, x, v )

Landau damping initialisation function

[Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)


"""
function landau( ϵ, kx, x, v )
    
    (1.0 .+ ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.*v))
    
end

#md # ---

nx, nv = 256, 256
xmin, xmax =  0., 4*π
vmin, vmax = -6., 6.
tf = 60
nt = 600
meshx = UniformMesh(xmin, xmax, nx)
meshv = UniformMesh(vmin, vmax, nv);

#md # --
            
# Initialize distribution function
x = meshx.points
v = meshv.points
ϵ, kx = 0.001, 0.5

# Allocate arrays for distribution function and its transposed
f = zeros(Complex{Float64},(nx,nv))
fᵀ= zeros(Complex{Float64},(nv,nx))
    
f .= landau( ϵ, kx, x, v); 

#md # ---

#nb Plot the distribution

#nb surface( x, v, real(f))

#nb #-
    
transpose!(fᵀ,f)
    
ρ  = compute_rho(meshv, f)
e  = zeros(ComplexF64, nx)
e .= compute_e(meshx, ρ)
    
mod_e = Float64[]
    
dt = tf / nt
    
advection_x! = AmpereAdvection( meshx )
advection_v! = AmpereAdvection( meshv );

#md # ---
        
@showprogress 1 for i in 1:nt
    
    advection_v!(fᵀ, e, 0.5dt)
    
    transpose!(f, fᵀ)
    
    advection_x!( f, e, v, dt)
    
    push!(mod_e, log(sqrt((sum(e.^2))*meshx.step))) # diagnostic
    
    transpose!(fᵀ, f)
    
    advection_v!(fᵀ, e, 0.5dt)

end

#md # ---

t = range(0, stop=tf, length=nt)
plot(t, -0.1533*t.-5.48)
plot!(t, mod_e , label=:ampere )
#md savefig("mod_e.svg")

@test true #src

#md # ![](mod_e.svg)

#md # ---
