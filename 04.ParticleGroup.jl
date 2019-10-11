# # Metaprogramming 
#
# ## The ParticleGroup example

using Distributions, Plots, Sobol, LinearAlgebra

#     ParticleGroup{D,V}(n_particles, charge, mass)
#
# - `D` : number of dimension in physical space
# - `V` : number of dimension in phase space
# - `n` : number of particles
# 
mutable struct ParticleGroup{D,V}

    n_particles       :: Int64
    data              :: Array{Float64, 2}

    function ParticleGroup{D,V}(n) where {D, V}

        data = zeros( Float64, (D+V, n))
        new( n, data) 
    end
end

# ---

# 
# Set position of ith particle of p to x
# 
@generated function set_x!( p :: ParticleGroup{D,V}, i, x :: Float64 ) where {D, V}

    :(p.data[1, i] = x)

end

# --

# 
# Set position of ith particle of p to x when x is a vector
# 
@generated function set_x!( p :: ParticleGroup{D,V}, i, x :: Vector{Float64} ) where {D, V}

     :(for j in 1:$D p.data[j, i] = x[j] end)

end

# ---

#
# Set velocity of ith particle of p to v
#
@generated function set_v!( p :: ParticleGroup{D,V}, i, v :: Float64 ) where {D, V}

    :(p.data[$D+1, i] = v)

end

# --

# 
# Set velocity of ith particle of p to v
# 
@generated function set_v!( p :: ParticleGroup{D,V}, i, v :: Vector{Float64} ) where {D, V}

    :(for j in 1:$V p.data[$D+j, i] = v[j] end)

end

# ---

#
# Get position of ith particle of p
#
@generated function get_x( p :: ParticleGroup{D,V}, i ) where {D, V}

     :(p.data[1:$D, i])

end

# 
# Get velocity of ith particle of p
# 
@generated function get_v( p :: ParticleGroup{D,V}, i ) where {D, V}

     :(p.data[$D+1:$D+$V, i])
end

# ---

# 
# Sampling from a probability distribution to initialize
# a Landau damping
# 
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
    
    s = SobolSeq(2)
    nbpart = pg.n_particles

    for i=1:nbpart

        v = sqrt(-2 * log( (i-0.5)/nbpart))
        r1, r2 = next!(s)
        θ = r1 * 2π
        set_x!(pg,  i, newton(r2))
        set_v!(pg,  i, [v * cos(θ), v * sin(θ)])
    end

end

# ---

n_particles = 10000
pg = ParticleGroup{1,2}( n_particles)
alpha, kx = 0.1, 0.5
landau_sampling!(pg, alpha, kx)

# --
xp = vcat([get_x(pg,i) for i in 1:pg.n_particles]...)
vp = vcat([get_v(pg, i) for i in 1:pg.n_particles]'...)

# ---

p = plot(layout=(3,1))
histogram!(p[1,1], xp, normalize=true, bins = 100, lab=:x)
plot!(p[1,1], x-> (1+alpha*cos(kx*x))/(2π/kx), 0., 2π/kx, lab="")
histogram!(p[2,1], vp[:,1], normalize=true, bins = 100, lab=:vx)
plot!(p[2,1], v-> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
histogram!(p[3,1], vp[:,2], normalize=true, bins = 100, lab=:vy)
plot!(p[3,1], v-> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
savefig("particles.svg")

# ![](particles.svg)

# ---

histogram2d(vp[:,1], vp[:,2], normalize=true, bins=100)
savefig("hist2d.svg")

# ![](hist2d.svg)

# ---
