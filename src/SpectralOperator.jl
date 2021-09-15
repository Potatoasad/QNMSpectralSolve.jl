abstract type SpectralOperator end

struct WaveEquationHyperboloidal{T1<:Function,T2<:Function,T3<:Function,T4<:Function} <: SpectralOperator
    w::T1
    p::T2
    qℓ::T3
    γ::T4
    x0::Float64
    x1::Float64
end

"""
    Defaults the compactified coordinate to be from -1 to 1
"""
function WaveEquationHyperboloidal(a::Function,b::Function,c::Function,d::Function; x0 = -1.0, x1 = 1.0)
    WaveEquationHyperboloidal(a,b,c,d,x0,x1)
end

ToChebyLine(t,x0,x1) = 2*(t-x0)/(x1-x0) - 1
ToOriginalLine(x,x0,x1) = 0.5*(x+1)*(x1-x0) + x0


"""
    Computes the diagonal Chebyshev matrix for a given function
    in position space
"""
function Diagonal(f::Function,N::Integer,x0,x1)
    xs = [cos(π*i/N) for i ∈ 0:N]
    diagonalvalues = f.(ToOriginalLine.(xs,x0,x1))
    diagm(diagonalvalues)
end


"""
    Computes the L1 Matrix of a particular wave equation in hyperboloidal coordinates
    acting on the basis of Chebyshev functions on the compactified coordinatep_func = L.p;
"""
function L1_matrix(L::WaveEquationHyperboloidal, N::Integer)
    #Save the functions to prevent repetitive accessing of functions from inside L
    w = L.w; p = L.p; qℓ = L.qℓ;

    #Transform functions to live on the standard Chebyshev [-1,1]
    x0 = L.x0; x1 = L.x1;
    p_cheby = ChebyshevT(p,Int64(N); x0=x0, x1=x1)
    ∂ₓp_cheby = ∂ₓ(p_cheby)
    ∂ₓp = (x -> ∂ₓp_cheby(x))

    #Define matrices that will be used to make the L1_matrix
    w⁻¹  =  Diagonal(x -> 1/w(x),N,x0,x1)
    pₘ   =  Diagonal(p,N,x0,x1);
    ∂ₓpₘ =  Diagonal(∂ₓp,N,x0,x1);
    qℓₘ  =  Diagonal(qℓ,N,x0,x1);
    ∂ₓₘ  =  (2/(x1-x0))*[DerivativeFunc(i,j,N) for i ∈ 0:N, j ∈ 0:N]

    #Construct the L1 matrix
    w⁻¹*(pₘ*∂ₓₘ*∂ₓₘ + ∂ₓpₘ*∂ₓₘ - qℓₘ)
end

"""
    Computes the L2 Matrix of a particular wave equation in hyperboloidal coordinates
    acting on the basis of Chebyshev functions on the compactified coordinatep_func = L.p;
"""
function L2_matrix(L::WaveEquationHyperboloidal, N::Integer)
    #Save the functions to prevent repetitive accessing of functions from inside L
    w = L.w; γ = L.γ;

    #Transform functions to live on the standard Chebyshev [-1,1]
    x0 = L.x0; x1 = L.x1;
    γ_cheby = ChebyshevT(γ,Int64(N); x0=x0, x1=x1)
    ∂ₓγ_cheby = ∂ₓ(γ_cheby)
    ∂ₓγ = (x -> ∂ₓγ_cheby(x))

    #Define matrices that will be used to make the L2_matrix
    w⁻¹  =  Diagonal(x -> 1/w(x),N,x0,x1)
    γₘ   =  Diagonal(γ,N,x0,x1);
    ∂ₓγₘ =  Diagonal(∂ₓγ,N,x0,x1);
    ∂ₓₘ  =  (2/(x1-x0))*[DerivativeFunc(i,j,N) for i ∈ 0:N, j ∈ 0:N]

    #println(w⁻¹)
    #println(γₘ)
    #println(∂ₓγₘ)
    #println(∂ₓₘ)
    #Construct the L2 matrix
    w⁻¹*(2*γₘ*∂ₓₘ + ∂ₓγₘ)
end

"""
    Computes the L operator (who's eigenvalues are the Quasinormal Modes)
"""
function matrix(F::SpectralOperator, N::Integer)
    L1 = L1_matrix(F,N);
    L2 = L2_matrix(F,N);
    -im*[zeros(N+1,N+1) I ; L1 L2]
end
