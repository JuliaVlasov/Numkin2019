# -*- coding: utf-8 -*-
# # Get more performance

# +
using Random, LinearAlgebra, BenchmarkTools

function test(A, B, C)
    C = C - A * B
    return C
end

A = rand(1024, 256)
B = rand(256, 1024)
C = rand(1024, 1024)
@btime test(A, B, C); #C, A and B are matrices. 
# -

function test_opt(A, B, C)
    BLAS.gemm!('N','N', -1., A, B, 1., C)
    return C
end
@btime test_opt(A, B, C); # avoids taking two unnecessary copies of the matrix C.

C = rand(1024, 1024)
all(test(A, B, C) .== test_opt(A, B, C))

# Memory alignement, and inplace computation.

# +
using FFTW

xmin, xmax, nx = 0, 4π, 1024
ymin, ymax, ny = 0, 4π, 1024
x = range(xmin, stop=xmax, length=nx+1)[1:end-1]
y = range(ymin, stop=ymax, length=ny+1)[1:end-1]
kx  = 2π ./ (xmax-xmin) .* [0:nx÷2-1;nx÷2-nx:-1]
ky  = 2π ./ (ymax-ymin) .* [0:ny÷2-1;ny÷2-ny:-1]
exky = exp.( 1im .* ky' .* x)

function df_dy( f, exky )
    ifft(exky .* fft(f, 2), 2)
end

f = sin.(x) .* cos.(y')

@btime df_dy(f, exky);

# +
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

@btime df_dy!(f, fᵗ, f̂ᵗ, exky );
# -

# # Pros
#
# - Packaging system is very efficient but i am not sure it will stay like this forever. The language is still young and when the number of package will grow...
# - PyPi 198,360 projects
# - R 14993 packages
# - Julia 3173 registered packages
# - It grows fast because it is very easy to create a package (easier than R and Python)
# - It is very easy to use GPU device.
# - Nice interface for Linear Algebra and Differential Equations
# - Easy access to BLAS and LAPACK

# # Cons
#
# - Julia is fast but it is not faster than Fortran. 
# - OpenMP is much better than the Julia multithreading library.
# - There is a MPI and PETSc package but they are not very active. 
# - The Julia community seems to prefer the ditributed processing approach. Like in Python community there has been an orientation towards data science in recent years. HPC is no longer in fashion and many of the latest developments are about machine learning and cloud computing.
# - Plotting takes time (20 seconds for the first plot)
# - Does not work well with vectorized code, you need to do a lot of inplace computation to avoid memory allocations and use explicit views to avoid copy.
# - Optimizing Julia code is often done at the expense of transparency. 

# # Julia is a language made for Science.
#
#  http://www.stochasticlifestyle.com/some-state-of-the-art-packages-in-julia-v1-0
#
#  * JuliaDiff – Differentiation tools
#  * JuliaDiffEq – Differential equation solving and analysis
#  * JuliaGeometry – Computational Geometry
#  * JuliaGraphs – Graph Theory and Implementation
#  * JuliaIntervals - Rigorous numerics with interval arithmetic & applications
#  * JuliaMath – Mathematics made easy in Julia
#  * JuliaOpt – Optimization
#  * JuliaPolyhedra – Polyhedral computation
#  * JuliaSparse – Sparse matrix solvers
#  * JuliaStats – Statistics and Machine Learning
#  * JuliaPlots - powerful convenience for visualization
#
