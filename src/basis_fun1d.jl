
abstract BasisFun
abstract BasisFun1d <: BasisFun
import Base.length
import Base.call



nmodes{T<:BasisFun}(b::T) = b.nmodes
isbndryint{T<:BasisFun}(b::T) = false
isbndryint{T<:BasisFun}(::Type{T}) = false

bndry_idx{T<:BasisFun}(b::T) = bndry_idx(b.lnum)
nbndry{T<:BasisFun}(b::T) = nbndry(b.lnum)

interior_idx{T<:BasisFun}(b::T) = bndry_idx(b.lnum)
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


immutable ModalC01d <: BasisFun1d
    nmodes::Int
    lnum::LocalNumSys
end
ModalC01d(n) = ModalC01d(n, LocalNumSys([1,2], [3:n;]))

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



immutable Lagrange1d{T<:Number} <: BasisFun1d
    nmodes::Int
    lnum::LocalNumSys
    z::Array{T,1}
end

function Lagrange1d{T<:Number}(n, ::Type{T}=Float64)
    z = Jacobi.zglj(n, 0, 0, T)
    Lagrange1d{T}(n, LocalNumSys([1,n], [2:n-1;]), z)
end

isbndryint(b::Lagrange1d) = true
basis1d{T<:Number}(b::Lagrange1d{T}, ξ::T, p) = Jacobi.lagrange(p, ξ, b.z)
