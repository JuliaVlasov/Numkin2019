---
layout: post
title:  Solving the Vlasov equation with Julia
---

*Numerical Methods for the Kinetic Equations of Plasma Physics*, Garching Oct 14-18, 2019.

The Vlasov equation is of fundamental importance in plasma physics
and particularly for simulations in Magnetic fusion and of Tokamak
plasmas. I present here a way to use the [Julia
language](https://julialang.org) to solve it numerically.

After a short introduction about the language, the first example
showed the kinetic simulation of Vlasov-Ampere system. The next
example proposes an implementation of the Particle In Cell method.

We are basing much of this effort on a previous implementation in
the Fortran language. We have found that the translation into Julia
is easy and it is interesting to look at what Julia has to offer
without degrading performance.

### Slides

1. [Introduction](/01/build/index.html) 
2. [Vlasov-Ampere with FFT](/02/build/index.html) 
3. [Rotation on GPU](/03/build/index.html) 
4. [Metaprogrammimg and Particle Group](/04/build/index.html) 
5. [GEMPIC](/05/build/index.html) 
6. [Conclusion](/06/build/index.html) 

### Jupyter notebooks

Either use the link above to open the notebooks in
[mybinder.org](https://mybinder.org/v2/gh/juliavlasov/Numkin2019/master?filepath=notebooks) or
run them locally:

```bash
git clone https://github.com/JuliaVlasov/Numkin2019
cd Numkin2019
julia --project
```

```julia
julia> using Pkg
julia> Pkg.instantiate()
julia> using IJulia
julia> notebook(dir=joinpath(pwd(),"notebooks"))
[ Info: running ...
```

