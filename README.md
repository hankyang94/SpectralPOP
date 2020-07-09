# SpectralPOP
- SpectralPOP is a Julia package of solving equality constrained polynomial optimization problems on a Euclidean sphere:

**inf_x { f(x) : x in R^n, hj(x) = 0, j = 1,...,l } with h1 := R - |x|^2,**

as well as extensive application in squared systems of polynomial equations solving:

**pj(x) = 0, j=1,...,n.**

- The main idea of SpectralPOP is to solve an SDP (moment) relaxation of the form:

**v = sup_X { <C,X> : X is psd, AX = b },**

which has constant trace property:

**AX = b => trace(X) = a,**

by using spectral (the largest eigenvalue) minimization:

**v = inf_z { aλ1(C - A'z) + b'z }**

with Limited memory bundle method instead of costly Interior-point methods.

- Compared to SumOfSquares (Mosek) and SketchyCGAL, SpectralPOP is cheaper, faster, but maintains the same accuracy with SumOfSquares on a tested sample of random dense equality constrained QCQPs on the unit sphere.

# Required softwares
SpectralPOP has been implemented on a desktop compute with the following softwares:
- Ubuntu 18.04.4
- Julia 1.3.1
- Fortran 2018

The following sofwares are used for comparison purposes:
- [SumOfSquares.jl](https://github.com/JuliaOpt/SumOfSquares.jl)
- [Mosek 9.1](https://www.mosek.com)
- [SketchyCGAL](https://github.com/alpyurtsever/SketchyCGAL)

# Remark
- Limited memory bundle method is only supported on Ubuntu.

# Installation
- To use SpectralPOP in Julia, run
```ruby
Pkg> add https://github.com/maihoanganh/SpectralPOP.git
```

# Usage
The following examples briefly guide to use SpectralPOP (see more examples in files test/test_....ipynb):

## Polynomial optimization
```ruby
using DynamicPolynomials

@polyvar x[1:2] # variables

f=x[1]^2+0.5*x[1]*x[2]-0.25*x[2]^2+0.75*x[1]-0.3*x[2] # objective function to minimize

R=1.0 # squared radius of a sphere constraint
h=[R-sum(x.^2);(x[1]-1.0)*x[2]] # equality constraints (including the sphere constraint)

k=2 # relaxed order

using SpectralPOP

# get the optimal value and an optimal solution
opt_val,opt_sol = CTP_POP(x,f,h,k,R,method="LMBM",EigAlg="Arpack",tol=1e-5)
```
Here:

- ```method="LMBM"```: [Limited memory bundle method](https://github.com/maihoanganh/LMBMinterface) (solver of spectral minimization). You can also choose ```method="PB"``` ([Proximal bundle method](https://github.com/maihoanganh/ProximalBundleMethod)) or ```method="SketchyCGAL"```.

- ```EigAlg="Arpack"```: [Arpack package](https://github.com/JuliaLinearAlgebra/Arpack.jl) (of computing the largest eigenvalue). You can also choose ```EigAlg="Normal"``` to use essential function of computing eigenvalue in Julia or ```method="Mix"``` to use the combination of the two methods.

- ```tol=1e-5```: the precision of the solver of spectral minimization.

See other options in the [link](https://github.com/maihoanganh/SpectralPOP/blob/master/examples/test_random_dense_quadratic_on_sphere.ipynb).

## Squared systems of polynomial equations

```ruby
using DynamicPolynomials

@polyvar x[1:2] # variables

# mickey equations
h=[x[1]^2 + 4*x[2]^2 - 4;
        2*x[2]^2 - x[1]]

L=10 # squared radius of a ball centered at origin containing at least one real root
k=1 # relaxed order

using SpectralPOP

# get a real root
root = ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-5)
```

# References
For more details, please refer to:

N. H. A. Mai, J.-B. Lasserre, and V. Magron. A hierarchy of spectral relaxations for polynomial optimization, 2020. Forthcoming.

(All numerical experiments of the paper are available in folder test/.)
