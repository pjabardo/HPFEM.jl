
"""
Quadrature information

Most spectral-element operations require calculation of integrals and derivatives.
The `QuadType` structure stores the required information when calculating integrals
and derivatives
"""
type QuadType{T<:Number}
    Q::Int
    z::Vector{T}
    w::Vector{T}
    D::Matrix{T}
end

"""
Constructs quadrature rules for Gauss-type quadrature.
"""
function QuadType{T<:Number,QT<:QUADRATURE_TYPE}(q::Int, ::Type{QT}=GLJ, ::Type{T}=Float64, a=0, b=0)
    z = Jacobi.qzeros(QT, q, a, b, T)
    w = Jacobi.qweights(QT, z, a, b)
    D = Jacobi.qdiff(QT, z, a, b)
    QuadType(q, z, w, D)
end

"""
Creates the ith Lagrange polynomial from a vector of nodes

Explicit construction of Lagrange polynomials from the nodes.
This function builds a `Poly` object.

 * `i` Index specifying which Lagrange polynomial should be returned.
 * `z` Vector containing the nodes of the Lagrange polynomials.
"""
function lagrange_poly{T<:Number}(i::Integer, z::AbstractVector{T})
    np = length(z)

    y = zeros(T, np-1)

    for j = 1:(i-1)
        y[j] = z[j]
    end
    for j = (i+1):np
        y[j-1] = z[j]
    end
    den = one(T)
    for j = 1:(np-1)
        den = den * (z[i] - y[j])
    end
    p = poly(y)
    return p/den
end

"""
Build Lagrange basis from any nodes

Probably not the fastest algorithm. But this won't be seriously used.
"""
function QuadType{T<:Number}(x::AbstractVector{T})  # Generic numbers:
    np = length(x)
    z = zeros(T,np)
    for i = 1:np
        z[i] = x[i]
    end

    w = zeros(T, np)
    D = zeros(T, np, np)

    for i = 1:np
        p = lagrange_poly(i, z)
        pint = polyint(p)
        pder = polyder(p)
        w[i] = polyval(pint, one(T)) - polyval(pint, -one(T))
        # Compute the rows of the derivative matrix
        for k = 1:np
            D[k,i] = polyval(pder, z[k])
        end
    end
    QuadType(np, z, w, D)
end
