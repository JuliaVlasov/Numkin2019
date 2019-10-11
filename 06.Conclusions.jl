# ### Optimizing Julia code is often done at the expense of transparency

using Random, LinearAlgebra, BenchmarkTools

function test(A, B, C)
    C = C - A * B
    return C
end

A = rand(1024, 256); B = rand(256, 1024); C = rand(1024, 1024)
@btime test(A, B, C); #C, A and B are matrices. 
nothing # hide

# --

function test_opt(A, B, C)
    BLAS.gemm!('N','N', -1., A, B, 1., C)
    return C
end
@btime test_opt(A, B, C); # avoids taking two unnecessary copies of the matrix C.
nothing # hide

# --

C = rand(1024, 1024)
all(test(A, B, C) .== test_opt(A, B, C))

# ---

# ### Derivative computation with FFT

using FFTW

xmin, xmax, nx = 0, 4π, 1024
ymin, ymax, ny = 0, 4π, 1024

x = range(xmin, stop=xmax, length=nx+1)[1:end-1]
y = range(ymin, stop=ymax, length=ny+1)[1:end-1]

ky  = 2π ./ (ymax-ymin) .* [0:ny÷2-1;ny÷2-ny:-1]
exky = exp.( 1im .* ky' .* x)

function df_dy( f, exky )
    ifft(exky .* fft(f, 2), 2)
end

f = sin.(x) .* cos.(y') # f is a 2d array created by broadcasting

@btime df_dy(f, exky);
nothing # hide

# ---

# ### Memory alignement, and inplace computation.

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

@btime df_dy!(f, fᵗ, f̂ᵗ, exky )
nothing # hide

# ---

# # Why use Julia language!
#
# - **You develop in the same language in which you optimize.**
# - Packaging system is very efficient (3173 registered packages)
# - PyPi (198,360 projects) R (14993 packages)
# - It is very easy to create a package (easier than R and Python)
# - It is very easy to use GPU device.
# - Nice interface for Linear Algebra and Differential Equations
# - Easy access to BLAS and LAPACK
# - Julia talks to all major Languages - mostly without overhead!

# ---

# # What's bad
#
# - It is still hard to build shared library or executable from Julia code.
# - Compilation times can be unbearable.
# - Plotting takes time (20 seconds for the first plot)
# - OpenMP is better than the Julia multithreading library but it is progressing.
# - There is a MPI and PETSc package but they are not very active. 
# - For parallelization, The Julia community seems to prefer the ditributed processing approach. 
# - Does not work well with vectorized code, you need to do a lot of inplace computation to avoid memory allocations and use explicit views to avoid copy.
# - Julia website proclaims that it is faster than Fortran but this is not true.

# [What's Bad About Julia by Jeff Bezanson](https://www.youtube.com/watch?v=TPuJsgyu87U)

# ---

# ## So when should i use Julia?
# 
# - Now! If you need performance and plan to write your own libraries.
# - In ~1-2 Years if you want a smooth deploy.
# - In ~3-5 Years if you want a 100% smooth experience.

# ## Julia Munich Meetup 
# Every two months, poll for the next meeting that
# will take place at Garching Campus : https://doodle.com/poll/z3ft2dytnaebyhh7.
#
# ## Python-Julia benchmarks by Thierry Dumont
#
# https://github.com/Thierry-Dumont/BenchmarksPythonJuliaAndCo/wiki

# ---

# # Julia is a language made for Science.
#
#  [Some State of the Art Packages](http://www.stochasticlifestyle.com/some-state-of-the-art-packages-in-julia-v1-0)
#
#  * JuliaDiffEq – Differential equation solving and analysis.
#  * JuliaDiff – Differentiation tools.
#  * JuliaGeometry – Computational Geometry.
#  * JuliaGraphs – Graph Theory and Implementation.
#  * JuliaIntervals - Rigorous numerics with interval arithmetic & applications.
#  * JuliaMath – Mathematics made easy in Julia.
#  * JuliaOpt – Optimization.
#  * JuliaPolyhedra – Polyhedral computation.
#  * JuliaSparse – Sparse matrix solvers.
#  * JuliaStats – Statistics and Machine Learning.
#  * JuliaPlots - powerful convenience for visualization.
#  * JuliaGPU - GPU Computing for Julia.
#  * FluxML - The Elegant Machine Learning Stack.
#
