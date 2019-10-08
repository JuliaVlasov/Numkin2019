# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     comment_magics: false
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

using Plots, Sobol

"""
Input r is a random number ``\\in [0,1]``

```math
    f(x) = 1 + \\alpha cos(k x)
```
on some domain ``[0, 2\\pi/k]``

Solve the equation ``P(x)-r=0`` with Newton’s method

```math
    x^{n+1} = x^n – (P(x)-(2\\pi r / k)/f(x) 
```

with 
```math
P(x) = \\int_0^x (1 + \\alpha cos(k y)) dy
```
```math
P(x) = x + \\frac{\\alpha}{k} sin (k x)
```
"""
function newton(r)
    x0, x1 = 0.0, 1.0
    alpha, k = 0.1, 0.5
    r *= 2 * pi / k
    while (abs(x1-x0) > 1e-12)
        p = x0 + alpha * sin( k * x0) / k 
        f = 1 + alpha * cos( k * x0)
        x0, x1 = x1, x0 - (p - r) / f
    end
    x1
end

?newton

# +
function landau( nbpart :: Int64)
    
   xp = Float64[]
   vp = Vector{Float64}[]
    
   s = SobolSeq(2)

   for k=1:nbpart

      v = sqrt(-2 * log( (k-0.5)/nbpart))
      r1, r2 = next!(s)
      θ = r1 * 2π
      push!(xp,  newton(r2))
      push!(vp,  [v * cos(θ), v * sin(θ)])

   end

   xp, vcat(vp'...)

end

xp, vp = landau(100000);
# -

p = plot(layout=(3,1))
histogram!(p[1,1], xp, normalize=true, bins = 100,  layout=(2,1), lab = "x")
plot!(p[1,1], x-> (1+0.1*cos(0.5*x))/4π, 0., 4π)
histogram!(p[2,1], vp[:,1], normalize=true, bins = 100,  layout=(2,1), lab = "vx")
plot!(p[2,1], x-> (exp(-x^2/2))/sqrt(2π), -6, 6)
histogram!(p[3,1], vp[:,2], normalize=true, bins = 100,  layout=(3,1), lab = "vy")
plot!(p[3,1], x-> (exp(-x^2/2))/sqrt(2π), -6, 6)
