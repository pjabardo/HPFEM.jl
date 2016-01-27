using Jacobi
using Base.LinAlg
abstract BasisFun
abstract BasisFun1d <: BasisFun
import Base.length
import Base.call

nmodes{T<:BasisFun1d}(b::T) = b.nmodes
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

import Base.convert
convert(::Type{LocalNumSys}, b::BasisFun) = b.lnum

type QuadType{T<:Number}
    Q::Int
    z::Vector{T}
    w::Vector{T}
    D::Matrix{T}
end
function QuadType{T<:Number,QT<:QUADRATURE_TYPE}(q::Int, ::Type{QT}=GLJ, ::Type{T}=Float64, a=0, b=0)
    z = qzeros(QT, q, a, b, T)
    w = qweights(QT, z, a, b)
    D = qdiff(QT, z, a, b)
    QuadType(q, z, w, D)
end



type Basis1d{T<:Number,B<:BasisFun1d}
    M::Int
    Q::Int
    ξ::Vector{T}
    w::Vector{T}
    D::Array{T,2}
    ϕ::Array{T,2}
    dϕ::Array{T,2}
    imass::Cholesky{T}
    bas::B
    quad::QuadType
end

function Basis1d{T<:Number, B<:BasisFun1d}(b::B, q::QuadType, ::Type{T}=Float64)
    # Create the basis function
    m = nmodes(b)
    # Obter as informações da quadratura
    Q = q.Q
    ξ = q.z
    w = q.w
    D = q.D
    ϕ = zeros(T, Q, m)
    # Preencher as funções de base:
    
    for k = 1:m
        for i=1:Q
            ϕ[i,k] = b(ξ[i], k)
        end
    end
    
    # Calcular as derivadas:
    dϕ = D * ϕ
    
    # calcular a matrix de massa
    mass = zeros(m,m)
    for k = 1:m
        for i = k:m
            mm = 0.0
            for j = 1:Q
                mm += ϕ[j,i] * ϕ[j,k] * w[j]
            end
            mass[i,k] = mm
            mass[k,i] = mm
        end
    end
    
    imass = cholfact(mass)
    Basis1d{T,B}(m, Q, ξ, w, D, ϕ, dϕ, imass, b, q)
end


qnodes(b::Basis1d) = b.ξ

nmodes(b::Basis1d) = b.M
nquad(b::Basis1d) = b.Q
basis_order(b::Basis1d) = b.M-1
weights(b::Basis1d) = b.w
basis(b::Basis1d) = b.ϕ
dbasis(b::Basis1d) = b.dϕ




basis1d(b::Basis1d, x, p) = basis1d(b.bas, x, p)

basis1d!(b::Basis1d, x::AbstractArray, y::AbstractArray, p) = basis1d!(b.bas, x, y, p)
basis1d(b::Basis1d, x::AbstractArray, p) = basis1d!(b, x, similar(x), p)

call(b::Basis1d, x, p) = basis1d(b.bas, x, p)

Basis1d(m::Int, q::Int) = Basis1d(ModalC01d(m), QuadType(q))
Basis1d(m::Int) = Basis1d(m, m+1)

nbndry(b::Basis1d) = nbndry(b.bas)
ninterior(b::Basis1d) = ninterior(b.bas)

bndry_idx(b::Basis1d) = bndry_idx(b.bas)
interior_idx(b::Basis1d) = interior_idx(b.bas)

seq2bi!{T}(b::Basis1d, x::AbstractVector{T}, y::AbstractVector{T}) = seq2bi!(b.bas.lnum, x, y)
seq2bi{T}(b::Basis1d, x::AbstractVector{T}) = seq2bi(b.bas.lnum, x)
bi2seq!{T}(b::Basis1d, x::AbstractVector{T}, y::AbstractVector{T}) = bi2seq!(b.bas.lnum, x, y)
bi2seq{T}(b::Basis1d, x::AbstractVector{T}) = bi2seq(b.bas.lnum, x)
bi2seq!{T}(b::Basis1d, xb::AbstractVector{T}, xi::AbstractVector{T}, y::AbstractVector{T}) =
    bi2seq!(b.bas.lnum, xb, xi, y)
bi2seq{T}(b::Basis1d, xb::AbstractVector{T}, xi::AbstractVector{T}) = bi2seq(b.bas.lnum, xb, xi)
seq2b!{T}(b::Basis1d, x::AbstractVector{T}, y::AbstractVector{T}) = seq2b!(b.bas.lnum, x, y)
seq2i!{T}(b::Basis1d, x::AbstractVector{T}, y::AbstractVector{T}) = seq2i!(b.bas.lnum, x, y)
seq2b{T}(b::Basis1d, x::AbstractVector{T}) = seq2b(b.bas.lnum, x)
seq2i{T}(b::Basis1d, x::AbstractVector{T}) = seq2i(b.bas.lnum, x)


function project(b::Basis1d, f::AbstractVector)

  ϕ = basis(b)
  w = weights(b)
  Q = nquad(b)
  M = nmodes(b)
  iM = b.imass
  fh = zeros(M)

  for k = 1:M
    F = 0.0
    for q = 1:Q
      F += f[q] * ϕ[q,k] * w[q]
    end
    fh[k] = F
  end

  A_ldiv_B!(iM, fh)

  return fh

end


project(b::Basis1d, f::Function) = project(b, f(qnodes(b)))


function mass_matrix(b::Basis1d)
  M = nmodes(b)
  Q = nquad(b)
  ϕ = basis(b)
  w = weights(b)
  mass = zeros(M,M)

  for k = 1:M
    for i = k:M
      m = 0.0
      for q = 1:Q
        m += ϕ[q,k] * ϕ[q, i] * w[q]
      end
      mass[k,i] = m
      mass[i,k] = m
    end
  end
  return mass
end




function stiff_matrix(b::Basis1d, mat)
  M = nmodes(b)
  Q = nquad(b)
  dϕ = dbasis(b)
  w = weights(w)

  mat = zeros(M,M)

  for k = 1:M
    for i = k:M
      m = 0.0
      for q = 1:Q
        m += dϕ[q,k] * dϕ[q, i] * w[q]
      end
      mat[k,i] = m
      mat[i,k] = m
    end
  end
  return mat
end















