using Jacobi

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

    new(m, q, ξ, w, D, ϕ, dϕ)
  end
end

num_modes(b::Basis1d) = b.M
num_quad(b::Basis1d) = b.Q
basis_order(b::Basis1d) = b.M-1


function modal_C0_basis(ξ, m, M)
    if m == 1
      ϕ = (1 - ξ) / 2
    elseif m == 2
      ϕ = (1 + ξ) / 2
    else
      ϕ = (1 - ξ)*(1 + ξ) / 4 * jacobi(ξ, m-2, 1, 1)
    end

    return ϕ
end


Basis1d(m, q) = Basis1d(m, q, modal_C0_basis)
Basis1d(m) = Basis1d(m, m+1)












