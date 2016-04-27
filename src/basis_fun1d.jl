
abstract BasisFun{T<:Number}
abstract BasisFun1d{T<:Number} <: BasisFun{T}

" Total number of modes of a basis function "
nmodes{T<:BasisFun}(b::T) = b.nmodes


" Total number of modes of a basis function "
length{T<:BasisFun1d}(b::T) = nmodes(b)

" Total number of modes of a basis function "
size{T<:BasisFun1d}(b::T) = (nmodes(b),)

" Can the modes be decoupled into boundary/interior modes? "
isbndryint{T<:BasisFun}(b::T) = false
isbndryint{T<:BasisFun}(::Type{T}) = false


bndidx{T<:BasisFun}(b::T) = bndidx(b.lnum)
nbndry{T<:BasisFun}(b::T) = nbndry(b.lnum)

intidx{T<:BasisFun}(b::T) = intidx(b.lnum)
ninterior{T<:BasisFun}(b::T) = ninterior(b.lnum)

seq2bi!{T}(b::BasisFun, x::AbstractVector{T}, y::AbstractVector{T}) = seq2bi!(b.lnum, x, y)
seq2bi{T}(b::BasisFun, x::AbstractVector{T}) = seq2bi(b.lnum, x)
bi2seq!{T}(b::BasisFun, x::AbstractVector{T}, y::AbstractVector{T}) = bi2seq!(b.lnum, x, y)
bi2seq{T}(b::BasisFun, x::AbstractVector{T}) = bi2seq(b.lnum, x)
bi2seq!{T}(b::BasisFun, xb::AbstractVector{T}, xi::AbstractVector{T}, y::AbstractVector{T}) =
    bi2seq!(b.lnum, xb, xi, y)
bi2seq{T}(b::BasisFun, xb::AbstractVector{T}, xi::AbstractVector{T}) = bi2seq(b.lnum, xb, xi)
seq2b!{T}(b::BasisFun, x::AbstractVector{T}, y::AbstractVector{T}) = seq2b!(b.lnum, x, y)
seq2i!{T}(b::BasisFun, x::AbstractVector{T}, y::AbstractVector{T}) = seq2i!(b.lnum, x, y)
seq2b{T}(b::BasisFun, x::AbstractVector{T}) = seq2b(b.lnum, x)
seq2i{T}(b::BasisFun, x::AbstractVector{T}) = seq2i(b.lnum, x)



"""
Compute the value of mode `i` at all points ξ

The values are returned on array `y`.
"""
function basis1d!{T<:BasisFun1d}(b::T, ξ::AbstractArray, y::AbstractArray, p::Integer)
    for i = 1:length(ξ)
        y[i] = basis1d(b, ξ[i], p)
    end
    y
end


basis1d{T<:BasisFun1d}(b::T, ξ::AbstractArray, p::Integer) = basis1d!(b, ξ, similar(ξ), p)
    
call{T<:BasisFun1d, N<:Number}(b::T, ξ::N, p::Integer) = basis1d(b, ξ, p)
call{T<:BasisFun1d, N<:Number}(b::T, ξ::AbstractArray{N}, p::Integer) = basis1d(b, ξ, p)

"""
Karniadakis/Sherwin 1D modal polynomial basis.
"""
immutable ModalC01d{T<:Number} <: BasisFun1d{T}
    nmodes::Int
    lnum::LocalNumSys1d
end
ModalC01d(n) = ModalC01d(n, LocalNumSys1d([1,2], [3:n;]))

isbndryint(b::ModalC01d) = true
isbndryint(::Type{ModalC01d}) = true

function basis1d{T<:Number}(b::ModalC01d, ξ::T, p::Integer)

    if p == 1
      ϕ = (one(T) - ξ) / 2
    elseif p == 2
      ϕ = (one(T) + ξ) / 2
    else
      ϕ = (one(T) - ξ)*(one(T) + ξ) / 4 * jacobi(ξ, p-3, 1, 1)
    end
    return ϕ
end

"""
Legendre orthogonal basis
"""
immutable Legendre1d{T<:Number} <: BasisFun1d{T}
    nmodes::Int
    lnum::LocalNumSys1d
end
Legendre1d(n) = Legendre1d(n+1, LocalNumSys1d( [1:(n+1);], Int[]) )

basis1d{T<:Number}(b::Legendre1d, ξ::T, p::Integer) = legendre(ξ, p-1)
    
"""
Lagrange polynomial basis.
"""
immutable Lagrange1d{T<:Number} <: BasisFun1d
    nmodes::Int
    lnum::LocalNumSys1d
    z::Array{T,1}
end

"""
Build a Lagrange polynomial basis at Gauss-Lobatto nodes.
"""
function Lagrange1d{T<:Number}(n, ::Type{T}=Float64)
    z = Jacobi.zglj(n, 0, 0, T)
    Lagrange1d{T}(n, LocalNumSys1d([1,n], [2:n-1;]), z)
end

"""
Build a Lagrange polynomial basis at any set of distinct nodes.
"""
function Lagrange1d{T<:Number}(x::AbstractVector{T})
    n = length(x)
    z = zeros(T, n)

    for i = 1:n
        z[i] = x[i]
    end
    bidx = Int[]
    iidx = Int[]
    
    for i = 1:n
        if z[i]==-1 || z[i]==1
            append!(bidx, [i])
        else
            append!(iidx, [i])
        end
    end
    Lagrange1d(n, LocalNumSys1d(bidx, iidx), z)
end


isbndryint(b::Lagrange1d) = isbndryint(b.lnum)

basis1d{T<:Number}(b::Lagrange1d{T}, ξ::T, p) = Jacobi.lagrange(p, ξ, b.z)
