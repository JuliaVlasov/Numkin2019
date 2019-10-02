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


# Pros

- Packaging system is very efficient but i am not sure it will stay like this forever. The language is still young and when the number of package will grow...
- PyPi 198,360 projects
- R 14993 packages
- Julia 3173 registered packages
- It grows fast because it is very easy to create a package (easier than R and Python)
- It is very easy to use GPU device.

# Cons
- Julia is fast but it is not faster than Fortran. OpenMP is much better than 
the multithreading 
- There is a MPI and PETSc package but they are not very active. The Julia community seems to prefer the ditributed processing approach they like to use in Data science.
