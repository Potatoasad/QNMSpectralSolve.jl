struct ChebyshevT{T}
    coeffs::Vector{T}
    N::Int64
    x0::Float64
    x1::Float64
end

ChebyshevT(c::AbstractArray;x0=-1.0,x1=1.0) = ChebyshevT(c,length(c)-1,x0,x1)

"""
    Functions to convert between compactified coordinates and the Chebyshev default [-1,1]
"""
ToChebyLine(t,x0,x1) = 2*(t-x0)/(x1-x0) - 1
ToOriginalLine(x,x0,x1) = 0.5*(x+1)*(x1-x0) + x0


"""
    ChebyshevT function evaluation
"""
function (f::ChebyshevT)(t)
    x = ToChebyLine(t,f.x0,f.x1)
    z = acos(x)
    return (f.coeffs[1]/2 + sum(f.coeffs[k+1]*cos(k*z) for k ∈ 1:(f.N)))
end


"""
    Kronecker Delta
"""
δ(i::Integer,j::Integer) = Int64(i == j)

"""
    kth Chebyshev Polynomial Evaluation
"""
T(k,x) = cos(k*acos(x))


"""
    Constructing a ChebyshevT from a given function g(x) over the domain [x0,x1]
"""
function ChebyshevT(g::Function, N::Integer; x0=-1.0, x1=1.0)
    x = [cos(π*i/N) for i ∈ 0:N]
    f(x) = g(ToOriginalLine(x,x0,x1))
    cᵢ = [(((2-δ(i,N))/(2*N))*(f(x[0+1])+((-1)^(i))*f(x[N+1])+
            2*sum((f(x[j+1])*T(i,x[j+1])) for j ∈ 1:(N-1)))) for i ∈ 0:N]
    ChebyshevT(cᵢ;x0=x0,x1=x1)
end


integral_coeffs(n_index) = n_index != 2 ? ((-1)^(n_index-1) + 1)/(1-(n_index-1)^2) : 0.0

"""
    Integrate a ChebyshevT function
"""
function Integrate(C::ChebyshevT)
    (C.coeffs[1]+sum(C.coeffs[k]*integral_coeffs(k) for k ∈ 2:(C.N+1)))*(C.x1-C.x0)/2
end



α(i,N) = (i==0 | i==N) ? 2 : 1

"""
    Matrix elements of the Derivative Operator on the Chebyshev grid
"""
function DerivativeFunc(i,j,N)
    if (i==N) & (j==N)
        return (-(2*N^2 + 1)/6)
    elseif (i==0) & (j==0)
        return ((2*N^2 + 1)/6)
    elseif (i>0)&(j<N)&(i==j)
        xj = cos(π*j/N)
        return -xj/(2*(1-xj^2))
    elseif (i≠j)
        xj = cos(π*j/N)
        xi = cos(π*i/N)
        return (α(i,N)/α(j,N))*((-1)^(abs(i-j))/(xi-xj))
    else
        return convert(typeof(N*0.1),0.0)
    end
end


"""
    Constructing a diagonal matrix from a function
"""
function MatrixfromFunc(f::Function, N::Integer)
    x = [cos(π*i/N) for i ∈ 0:N]
    diagm(f.(x))
end

"""
    Chebyshev Polynomials of the 2nd kind. Mostly for the intermediate representation of the derivative
    of the Chebyshev polynomials of the 1st kind. These are then interpolated on the grid and converted
    back to polynomials of the first kind
"""
struct ChebyshevU{T}
    coeffs::Vector{T}
    N::Int64
    x0::Float64
    x1::Float64
end

"""
    Evaluate a ChebyshevU interpolant
"""
function (f::ChebyshevU)(t)
    x = ToChebyLine(t,f.x0,f.x1)
    if x == -1.0
        return (f.coeffs[1]/2 + sum(f.coeffs[k+1]*(k+1)*(-1)^k for k ∈ 1:(f.N)))
    elseif x == 1.0
        return (f.coeffs[1]/2 + sum(f.coeffs[k+1]*(k+1) for k ∈ 1:(f.N)))
    else
        θ = acos(x)
        return (f.coeffs[1]/2 + sum(f.coeffs[k+1]*sin((k+1)*θ) for k ∈ 1:(f.N))/(sin(θ)))
    end
end


"""
    Calculate the derivative of a ChebyshevT polynomial by first going through the intermediate representation
    of a ChebyshevU polynomial.
    Step 1. Convert coefficients to ChebyshevU coefficients
    Step 2. Sample the ChebyshevU interpolant on the ChebyshevT grid
    Step 3. Use the samples to get the new ChebyshevT interpolant
"""
function ∂ₓ(f::ChebyshevT)
    newcoeffs = (2/(f.x1-f.x0))*[2*f.coeffs[2] ; [f.coeffs[n+2]*(n+1) for n ∈ 1:(f.N-1)]... ; zero(f.coeffs[end])]
    fU = ChebyshevU(newcoeffs,f.N,f.x0,f.x1)
    ChebyshevT(x->fU(x), f.N; x0=f.x0, x1=f.x1)
end
