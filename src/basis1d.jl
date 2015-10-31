using Jacobi
using Base.LinAlg
abstract BasisFun
abstract BasisFun1d <: BasisFun
import Base.length
import Base.call

nmodes{T<:BasisFun1d}(b::T) = b.nmodes
isbndryint{T<:BasisFun}(b::T) = false
isbndryint{T<:BasisFun}(::Type{T}) = false

bndry_idx{T<:BasisFun}(b::T) = b.bndry
nbndry{T<:BasisFun}(b::T) = length(b.bndry)

interior_idx{T<:BasisFun}(b::T) = b.interior
ninterior{T<:BasisFun}(b::T) = length(b.interior)

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
    bndry::Array{Int,1}
    interior::Array{Int,1}
end
ModalC01d(n) = ModalC01d(n, [1,2], [3:n;])

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
    bndry::Array{Int,1}
    interior::Array{Int,1}
    z::Array{T,1}
end

function Lagrange1d{T<:Number}(n, ::Type{T}=Float64)
    z = Jacobi.zglj(n, 0, 0, T)
    Lagrange1d{T}(n, [1,n], [2:n-1;], z)
end

isbndryint(b::Lagrange1d) = true
basis1d{T<:Number}(b::Lagrange1d{T}, ξ::T, p) = Jacobi.lagrange(p, ξ, b.z)
    
type Basis1d{T<:Number,B<:BasisFun1d}
    M::Int
    Q::Int
    ξ::Vector{T}
    w::Vector{T}
    D::Array{T,2}
    ϕ::Array{T,2}
    dϕ::Array{T,2}
    imass::Cholesky{T}
    b::B
end

function Basis1d{T<:Number, B<:BasisFun1d}(b::B, q::Integer, ::Type{T}=Float64)
        # Create the basis function
        m = nmodes(b)
        # Obter as informações da quadratura
        ξ = zglj(q)
        w = wglj(ξ)
        D = dglj(ξ)
        ϕ = zeros(q, m)
        # Preencher as funções de base:
        
        for k = 1:m
            for i=1:q
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
                for j = 1:q
                    mm += ϕ[j,i] * ϕ[j,k] * w[j]
                end
                mass[i,k] = mm
                mass[k,i] = mm
            end
        end
        
        imass = cholfact(mass)
        Basis1d{T,B}(m, q, ξ, w, D, ϕ, dϕ, imass, b)
end


qnodes(b::Basis1d) = b.ξ

num_modes(b::Basis1d) = b.M
num_quad(b::Basis1d) = b.Q
basis_order(b::Basis1d) = b.M-1
weights(b::Basis1d) = b.w
basis(b::Basis1d) = b.ϕ
dbasis(b::Basis1d) = b.dϕ




basis1d(b::Basis1d, x, p) = basis1d(b.b, x, p)

basis1d!(b::Basis1d, x::AbstractArray, y::AbstractArray, p) = basis1d!(b.b, x, y, p)
basis1d(b::Basis1d, x::AbstractArray, p) = basis1d!(b, x, similar(x), p)

call(b::Basis1d, x, p) = basis1d(b.b, x, p)

Basis1d(m, q) = Basis1d(ModalC01d(m), q)
Basis1d(m) = Basis1d(m, m+1)


function project(b::Basis1d, f::AbstractVector)

  ϕ = basis(b)
  w = weights(b)
  Q = num_quad(b)
  M = num_modes(b)
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
  M = num_modes(b)
  Q = num_quad(b)
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
  M = num_modes(b)
  Q = num_quad(b)
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















