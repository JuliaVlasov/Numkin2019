ENV["GKSwstype"]="100" #src

# ### Optimizing Julia code is often done at the expense of transparency

using Random, LinearAlgebra, BenchmarkTools

function test(A, B, C)
    C = C - A * B
    return C
end

A = rand(1024, 256); B = rand(256, 1024); C = rand(1024, 1024)
@btime test(A, B, C); #C, A and B are matrices. 
#md nothing # hide

#md # --

function test_opt(A, B, C)
    BLAS.gemm!('N','N', -1., A, B, 1., C)
    return C
end
@btime test_opt(A, B, C) # avoids taking two unnecessary copies of the matrix C.
#md nothing # hide

#md # --

C = rand(1024, 1024)
all(test(A, B, C) .== test_opt(A, B, C))

#md # ---

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

@btime df_dy(f, exky)
#md nothing # hide

#md # ---

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
#md nothing # hide

#md # ---

#md # # Why use Julia language!
#md #
#md # - **You develop in the same language in which you optimize.**
#md # - Packaging system is very efficient (3173 registered packages)
#md # - PyPi (198,360 projects) R (14993 packages)
#md # - It is very easy to create a package (easier than R and Python)
#md # - It is very easy to use GPU device.
#md # - Nice interface for Linear Algebra and Differential Equations
#md # - Easy access to BLAS and LAPACK
#md # - Julia talks to all major Languages - mostly without overhead!

#md # ---

#md # # What's bad
#md #
#md # - It is still hard to build shared library or executable from Julia code.
#md # - Compilation times can be unbearable.
#md # - Plotting takes time (20 seconds for the first plot)
#md # - OpenMP is better than the Julia multithreading library but it is progressing.
#md # - There is a MPI and PETSc package but they are not very active. 
#md # - For parallelization, The Julia community seems to prefer the ditributed processing approach. 
#md # - Does not work well with vectorized code, you need to do a lot of inplace computation to avoid memory allocations and use explicit views to avoid copy.
#md # - Julia website proclaims that it is faster than Fortran but this is not true. But it is very close and it is progressing.
#md 
#md # [What's Bad About Julia by Jeff Bezanson](https://www.youtube.com/watch?v=TPuJsgyu87U)

#md # ---

#md # ## So when should i use Julia?
#md # 
#md # - Now! If you need performance and plan to write your own libraries.
#md # - In ~1-2 Years if you want a smooth deploy.
#md # - In ~3-5 Years if you want a 100% smooth experience.
#md 
#md # ## Julia Munich Meetup 
#md # Every two months, poll for the next meeting that
#md # will take place at Garching Campus : https://doodle.com/poll/z3ft2dytnaebyhh7.
#md #
#md # ## Python-Julia benchmarks by Thierry Dumont
#md #
#md # https://github.com/Thierry-Dumont/BenchmarksPythonJuliaAndCo/wiki
#md 
#md # ---
#md 
#md # # Julia is a language made for Science.
#md #
#md #  [Some State of the Art Packages](http://www.stochasticlifestyle.com/some-state-of-the-art-packages-in-julia-v1-0)
#md #
#md #  * JuliaDiffEq – Differential equation solving and analysis.
#md #  * JuliaDiff – Differentiation tools.
#md #  * JuliaGeometry – Computational Geometry.
#md #  * JuliaGraphs – Graph Theory and Implementation.
#md #  * JuliaIntervals - Rigorous numerics with interval arithmetic & applications.
#md #  * JuliaMath – Mathematics made easy in Julia.
#md #  * JuliaOpt – Optimization.
#md #  * JuliaPolyhedra – Polyhedral computation.
#md #  * JuliaSparse – Sparse matrix solvers.
#md #  * JuliaStats – Statistics and Machine Learning.
#md #  * JuliaPlots - powerful convenience for visualization.
#md #  * JuliaGPU - GPU Computing for Julia.
#md #  * FluxML - The Elegant Machine Learning Stack.
#md #

@test true #src
