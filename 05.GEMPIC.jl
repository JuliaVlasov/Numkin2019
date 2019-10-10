# # GEMPIC
#
# Geometric ElectroMagnetic Particle-In-Cell Methods
#
# https://arxiv.org/abs/1609.03053
#
# Michael Kraus, Katharina Kormann, Philip J. Morrison, Eric Sonnendrücker
#
# Framework for Finite Element Particle-in-Cell methods based on 
# the discretization of the underlying Hamiltonian structure of the 
# Vlasov-Maxwell system. 

# ```julia
# using Pkg
# Pkg.add(PackageSpec(url="https://github.com/juliavlasov/GEMPIC.jl"))
# ```

using ProgressMeter, Plots, GEMPIC

#---

# # Strong Landau Damping
#
# The physical parameters 

kx, α = 0.5, 0.5
xmin, xmax = 0, 2π/kx
domain = [xmin, xmax, xmax - xmin]

# --

# The numerical parameters

∆t = 0.05
nx = 32 
n_particles = 100000
mesh = Mesh( xmin, xmax, nx)
spline_degree = 3

# --

# Initialize particles

mass, charge = 1.0, 1.0
particle_group = ParticleGroup{1,2}( n_particles, mass, charge, 1)   
sampler = LandauDamping( α, kx)

sample!(sampler, particle_group)

# ---

# Particle-mesh coupling operators

kernel_smoother1 = ParticleMeshCoupling( domain, [nx], n_particles, 
                                         spline_degree-1, :galerkin)    

kernel_smoother0 = ParticleMeshCoupling( domain, [nx], n_particles, 
                                         spline_degree, :galerkin)
# Allocate electrostatic fields

rho = zeros(Float64, nx)
ex  = zeros(Float64, nx);

# Maxwell solver

maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)

# Set initial electric field from ρ

solve_poisson!( ex, particle_group, 
                kernel_smoother0, maxwell_solver, rho)

# ---

# ### Charge density ρ(x)

xg = LinRange(xmin, xmax, nx)
sval = eval_uniform_periodic_spline_curve(spline_degree-1, rho)
plot( xg, sval, label="ρ")
savefig("rho.svg")

# ![](rho.svg)

# ---

# ### Electric field e(x)

efield_poisson = zeros(Float64, nx)
# Init!ialize the field solver
maxwell_solver = Maxwell1DFEM(domain, nx, spline_degree)
# efield by Poisson
solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, rho )
sval = eval_uniform_periodic_spline_curve(spline_degree-1, efield_poisson)
plot( xg, sval )       
savefig("ex.svg")

# ![](ex.svg)

# ---

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

# ---

# ## Loop over time

steps, Δt = 100, 0.05

@showprogress 1 for j = 1:steps # loop over time

    ## Strang splitting
    strang_splitting!(propagator, Δt, 1)

    ## Diagnostics
    solve_poisson!( efield_poisson, particle_group, 
                    kernel_smoother0, maxwell_solver, rho)
    
    write_step!(thdiag, j * Δt, spline_degree, 
                    efield_dofs,  bfield_dofs,
                    efield_dofs_n, efield_poisson)

end

# ---

# ## Diagnostics stored in a dataframe

using DataFrames
first(thdiag.data, 5)
#-

plot(thdiag.data[!,:Time], log.(thdiag.data[!,:PotentialEnergyE1]))
savefig("thdiag.svg")

# ![](thdiag.svg)
