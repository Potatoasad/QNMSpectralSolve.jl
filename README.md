# QNMSpectralSolve.jl
This package follows what is done in https://arxiv.org/abs/2004.06434 and sets up some wave equations in hyperboloidal slicing. 
It uses Spectral Chebyshev methods as outlined in the paper to cast the problem of getting Quasinormal Modes of a transformed wave equation to a problem of finding the eigenvalues of a matrix. 

The eigenvalue problem is not very numerically stable, hence one needs to use `BigFloat`s and use a large of Chebyshev basis of `~200`. 

```julia
using QNMSpectralSolve
using GenericLinearAlgebra
using Plots

## Following Eq 70. in arXiv:2004.06434v4
w(σ) = 2*(1+σ)
p(σ) = 2*(σ^2)*(1-σ)
ℓ = 2; s = 2;
Vℓ(σ) = 2*(ℓ*(ℓ+1) + (1-s^2)*σ)
γ(σ) = 1-2*σ^2

## Define the equation using the functions above and get the matrix
□ = WaveEquationHyperboloidal(w,p,Vℓ,γ; x0=0.0,x1=1.0)
L = matrix(□,big(40)) #Get bigfloats by plugging in big(N)

## Solve the eigenvalue problem
Es = eigvals(L2)
x = real.(Es)
y = imag.(Es)
scatter(x,y) |> display

```
