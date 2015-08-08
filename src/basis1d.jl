using Jacobi
using Base.LinAlg
abstract BasisFun
abstract BasisFun1d
import Base.length



type Basis1d
  M::Int
  Q::Int
  ξ::Vector{Float64}
  w::Vector{Float64}
  D::Array{Float64,2}
  ϕ::Array{Float64,2}
  dϕ::Array{Float64,2}
  imass::Cholesky{Float64}
  fun::Function

  function Basis1d(m, q, fun)
    # Obter as informações da quadratura
    ξ = zglj(q)
    w = wglj(ξ)
    D = dglj(ξ)
    ϕ = zeros(q, m)
    # Preencher as funções de base:
    for k = 1:m
      for i=1:q
        ϕ[i,k] = fun(ξ[i], k, m)
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
    new(m, q, ξ, w, D, ϕ, dϕ, imass, fun)
  end
end


qnodes(b::Basis1d) = b.ξ

num_modes(b::Basis1d) = b.M
num_quad(b::Basis1d) = b.Q
basis_order(b::Basis1d) = b.M-1
weights(b::Basis1d) = b.w
basis(b::Basis1d) = b.ϕ
dbasis(b::Basis1d) = b.dϕ




basis1d(b::Basis1d, x, m) = b.fun(x, m, b.M)
function basis1d!(b::Basis1d, x::AbstractArray, m, y::AbstractArray)
  for i = 1:length(x)
    y[i] = basis1d(b, x[i], m)
  end
  return y
end
basis1d(b::Basis1d, x::AbstractArray, m) = basis1d!(b, x, m, zeros(x))


function modal_C0_basis(ξ, m, M)
    if m == 1
      ϕ = (1 - ξ) / 2
    elseif m == 2
      ϕ = (1 + ξ) / 2
    else
      ϕ = (1 - ξ)*(1 + ξ) / 4 * jacobi(ξ, m-3, 1, 1)
    end

    return ϕ
end


Basis1d(m, q) = Basis1d(m, q, modal_C0_basis)
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




type Element1d
  basis::Basis1d
  a::Float64
  b::Float64
  J::Float64
  Element1d(bas::Basis1d, a=-1.0, b=1.0) = new(bas, a, b, 2/(b-a))
end

num_modes(el::Element1d) = el.basis.M
num_quad(el::Element1d) = el.basis.Q
basis1d(el::Element1d) = el.basis
basis(el::Element1d) = basis(basis1d(el))
dbasis(el::Element1d) = dbasis(basis1d(el))
weights(el::Element1d) = weights(basis1d(el))
jacobian(el::Element1d) = el.J

function qnodes(el::Element1d)
  ξ = qnodes(basis1d(el))
  (1 - ξ)*el.a/2 + (1 + ξ)*el.b/2
end


mass_matrix(el::Element1d) = jacobian(el) * mass_matrix(basis1d(el))
stiff_matrix(el::Element1d) = jacobian(el)^3 * stiff_matrix(basis1d(el))



project(el::Element1d, f::AbstractVector) = project(basis1d(el), f)
project(el::Element1d, f::Function) = project(basis1d(el), qnodes(el))











