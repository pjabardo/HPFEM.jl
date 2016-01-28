import Base.length
import Base.call








import Base.convert
convert(::Type{LocalNumSys}, b::BasisFun) = b.lnum


abstract GenBasis
abstract GenBasis1d

type Basis1d{T<:Number,B<:BasisFun1d} <: GenBasis1d
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















