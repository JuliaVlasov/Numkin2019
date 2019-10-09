# -*- coding: utf-8 -*-
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

# +
using Random, LinearAlgebra, BenchmarkTools

function test(A, B, C)
    C = C - A*B
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
# -

# Memory alignement, and inplace computation.

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

# # Introduction
#
# High-level languages like Python and R let one explore and experiment rapidly, but can run slow. They also provide a large choice of packages and it is easy to share your work without care of dependencies. With thess langages, end-users have 
# a better expenrience et it is easy to provide a nice documentation with
# your software. 
#
# Low-level languages like Fortran/C++
# tend to take longer to develop, but run fast. It is often hard to manage dependencies, testing, documentation and packaging if you want to share your work. Many good libraries tend to wrap C++/Fortran code with Python or R to facilitate the use and offer a better interface or documentation. This is sometimes
# called the "two language problem" and is something the Julia
# developers set out to eliminate.
#
# When you mix two langages you have to use building tools like
# Makefiles. It is not easy to configure, sometimes it is not portable.
#
# You can have bootlenecks because you need to convert types from one
# language to the other.  You could need more memory. In R and Python
# you have to vectorize your implementation.
#
# - Julia's promise is to provide a "best of both worlds" experience for 
# programmers who need to develop novel algorithms and bring them into 
# production environments with minimal effort.
# You develop in the same language in which you optimize.
