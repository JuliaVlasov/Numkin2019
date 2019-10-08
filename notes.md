# Who am I ?

 - My name is *Pierre Navaro*
 - **Fortran 77 + PVM** : during my PhD 1998-2002 (Université du Havre) 
 - **Fortran 90-2003 + OpenMP-MPI** : Engineer in Strasbourg (2003-2015) 
 - **Numpy + Cython, R + Rcpp** : engineer in Rennes (2015-now)
 - **Julia v1.0** since July 2018 

 # Advection equation for a rotation in two dimensional domain

 $$
 \frac{d f}{dt} +  (y \frac{d f}{dx} - x \frac{d f}{dy}) = 0
 $$

 $$ 
 x \in [-\pi, \pi],\qquad y \in [-\pi, \pi] \qquad \mbox{ and } \qquad t \in [0, 200\pi] 
 $$

# ---

# Introduction

High-level languages like Python and R let one explore and experiment rapidly, but can run slow. They also provide a large choice of packages and it is easy to share your work without care of dependencies. With thess langages, end-users have 
a better expenrience et it is easy to provide a nice documentation with
your software. 

Low-level languages like Fortran/C++
tend to take longer to develop, but run fast. It is often hard to manage dependencies, testing, documentation and packaging if you want to share your work. Many good libraries tend to wrap C++/Fortran code with Python or R to facilitate the use and offer a better interface or documentation. This is sometimes
called the "two language problem" and is something the Julia
developers set out to eliminate.

When you mix two langages you have to use building tools like
Makefiles. It is not easy to configure, sometimes it is not portable.

You can have bootlenecks because you need to convert types from one
language to the other.  You could need more memory. In R and Python
you have to vectorize your implementation.

- Julia's promise is to provide a "best of both worlds" experience for 
programmers who need to develop novel algorithms and bring them into 
production environments with minimal effort.
You develop in the same language in which you optimize.


# Pros

- Packaging system is very efficient but i am not sure it will stay like this forever. The language is still young and when the number of package will grow...
- PyPi 198,360 projects
- R 14993 packages
- Julia 3173 registered packages
- It grows fast because it is very easy to create a package (easier than R and Python)
- It is very easy to use GPU device.
- Nice interface for Linear Algebra and Differential Equations
- Easy access to BLAS and LAPACK

# Cons
- Julia is fast but it is not faster than Fortran. OpenMP is much better than 
the multithreading 
- There is a MPI and PETSc package but they are not very active. The Julia community seems to prefer the ditributed processing approach. 
Like in Python community There has been an orientation towards data science in recent years. HPC is no longer in fashion and now all the latest developments are about machine learning and cloud computing.
- Plotting takes time (20 seconds for the first plot)
- Does not work well with vectorized code, you need to do a lot of inplace computation
to avoid memory allocations and use explicit views to avoid copy.
- Optimizing Julia code is often done at the expense of transparency, it is obvious what 

 C = C - A*B

means when C, A and B are matrices.  It is less obvious what

  BLAS.gemm!('N','N',-1.,A,B,1.,C)

means but the second form avoids taking two unnecessary copies of the matrix C.

exky = exp.( 1im * dt .* ky' .* x)
f .= ifft(exky .* fft(f, 2), 2)

fᵗ = zeros(ComplexF64, reverse(size(f)))
f̂ᵗ = zeros(ComplexF64, reverse(size(f)))
FFTW.set_num_threads(4)
fft_plan = plan_fft(fᵗ, 1, flags=FFTW.PATIENT)

transpose!(fᵗ,f)
mul!(f̂ᵗ,  fft_plan, fᵗ)
f̂ᵗ .= f̂ᵗ .* exky
ldiv!(fᵗ, fft_plan, f̂ᵗ)
transpose!(f, fᵗ)

We save some time because of the memory alignement, and inplace computation.
