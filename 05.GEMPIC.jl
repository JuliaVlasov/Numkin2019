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

# ## GEMPIC
#
# Geometric ElectroMagnetic Particle-In-Cell Methods
#
# https://arxiv.org/abs/1609.03053

using Pkg
Pkg.add(PackageSpec(url="https://github.com/juliavlasov/GEMPIC.jl"))

using ProgressMeter, Plots, GEMPIC

# # Strong Landau Damping
#
# Electrostatic example of strong Landau damping
# $$
# f(x,v) =\frac{1}{2\pi\sigma^2}  \exp 
# \Big( - \frac{v_1^2 + v_2^2}{2\sigma^2} \Big)
# ( 1+\alpha \cos(k x)􏰁),
# $$
#
# The physical parameters 

# +
kx, α = 0.5, 0.5
xmin, xmax = 0, 2π/kx
domain = [xmin, xmax, xmax - xmin]
∆t = 0.05
nx = 32 
n_particles = 100000
mesh = Mesh( xmin, xmax, nx)
spline_degree = 3

mass, charge = 1.0, 1.0
particle_group = ParticleGroup{1,2}( n_particles, mass, charge, 1)   
sampler = LandauDamping( α, kx)
# -

sample!(sampler, particle_group)

# +
xp = vcat([get_x(particle_group, i) for i in 1:n_particles]...)
vp = vcat([get_v(particle_group, i) for i in 1:n_particles]'...)

p = histogram([xp, vp[:,1], vp[:,2]], normalize=true, bins = 100,  layout=(3,1), lab = "")
plot!(p[1,1], x-> (1+α*cos(kx*x))/(2π/kx), 0., 2π/kx, lab="")
plot!(p[2,1], v-> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
plot!(p[3,1], v-> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")

# +
kernel_smoother1 = ParticleMeshCoupling( domain, [nx], n_particles, spline_degree-1, :galerkin)    
kernel_smoother0 = ParticleMeshCoupling( domain, [nx], n_particles, spline_degree, :galerkin)
rho = zeros(Float64, nx)
ex = zeros(Float64, nx)
maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)
solve_poisson!( ex, particle_group, 
                kernel_smoother0, maxwell_solver, rho)

xg = LinRange(xmin, xmax, nx)
sval = eval_uniform_periodic_spline_curve(spline_degree-1, rho)
plot( xg, sval, label="ρ")
# -

efield_poisson = zeros(Float64, nx)
# Init!ialize the field solver
maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)
# efield by Poisson
solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )
sval = eval_uniform_periodic_spline_curve(spline_degree-1, efield_poisson)
plot( xg, sval )       

# +
# # +
# Initialize the arrays for the spline coefficients of the fields
efield_dofs = [efield_poisson, zeros(Float64, nx)]
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

steps, Δt = 500, 0.05

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
# -

using DataFrames
first(thdiag.data, 10)

plot(thdiag.data[!,:Time], log.(thdiag.data[!,:PotentialEnergyE1]))