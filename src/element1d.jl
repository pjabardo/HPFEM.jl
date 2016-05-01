
abstract Element


type Element1d{T <: Number} <: Element
    id::Int
    a::T
    b::T
    ξ::Vector{T}
    w::Vector{T}
    D::Matrix{T}
    x::Vector{T}
    J::Vector{T}
    wJ::Vector{T}
    dξdx::Vector{T}
end


function Element1d{T<:Number}(id, a::T, b::T, bas::GenBasis1d)
    w = qweights(bas)
    ξ = qnodes(bas)
    Q = nquad(bas)
    x = zeros(T, Q)
    dξdx = zeros(T,Q)
    J = zeros(T,Q)
    wJ = zeros(T,Q)
    D = diffmat(bas)
    if !isinf(a) && !isinf(b)
        d = (b-a) / 2
        
        for i = 1:Q
            J[i] = d
            x[i] = (one(T)-ξ[i])*a/2 + (one(T) + ξ[i])*b/2
            dξdx[i] = one(T)/d
            wJ[i] = w[i] * d
        end
    end

    Element1d(id, a, b, ξ, w, D, x, J, wJ, dξdx)
end


eid(e::Element) = e.id
jacweights(e::Element) = e.wJ
deriv_ξ(e::Element1d) = e.dξdx
jacobian(el) = el.J
diffmat(e::Element1d) = e.D




project(el::Element1d, f::AbstractVector) = project(basis1d(el), f)
project(el::Element1d, f::Function) = project(basis1d(el), qnodes(el))
