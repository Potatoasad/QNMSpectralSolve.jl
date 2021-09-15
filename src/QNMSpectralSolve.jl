module QNMSpectralSolve

using LinearAlgebra

include("Chebyshev.jl")
include("SpectralOperator.jl")


export WaveEquationHyperboloidal, matrix, L1_matrix, L2_matrix
export ChebyshevT, ∂ₓ, Integrate, SpectralOperator

end # module
