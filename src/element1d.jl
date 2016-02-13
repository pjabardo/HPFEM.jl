
abstract Element


type Element1d{T <: Number} <: Element
    id::Int
    a::T
    b::T
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
    
    if !isinf(a) && !isinf(b)
        d = (b-a) / 2
        
        for i = 1:Q
            J[i] = d
            x[i] = (1-ξ[i])*a/2 + (1 + ξ[i])*b/2
            dξdx[i] = 1/d
            wJ[i] = w[i] * d
        end
    end

    Element1d(id, a, b, x, J, wJ, dξdx)
end


eid(e::Element) = e.id
jweights(e::Element) = e.wJ
deriv_ξ(e::Element1d) = e.dξdx


function mass_matrix!(el::Elemen1d, bas::GenBasis1d, λ::GenFunction1d, A::AbstractMatrix)
    
    

end

function add_mass_matrix!{T}(el::Element1d, mass::AbstractMatrix{T}, λ=one(T))
    J = jacobian(el)
    M = nmodes(el)
    Q = nquad(el)
    ϕ = basis1d(el)
    w = weights(el)
    for k = 1:M
        for i = k:M
            m = 0.0
            for q = 1:Q
                m += ϕ[q,k] * ϕ[q, i] * w[q]
            end
            mass[k,i] += m*J*λ
            if k != i
                mass[i,k] += m*J*λ
            end
        end
    end
    return mass
end


function add_mass_matrix!{T}(el::Element1d, mass::AbstractMatrix{T}, λ::AbstractVector{T})
    J = jacobian(el)
    M = nmodes(el)
    Q = nquad(el)
    ϕ = basis1d(el)
    w = weights(el)
    for k = 1:M
        for i = k:M
            m = 0.0
            for q = 1:Q
                m += ϕ[q,k] * ϕ[q, i] * w[q] * λ[q]
            end
            mass[k,i] += m*J*λ
            if k != i
                mass[i,k] += m*J*λ
            end
        end
    end
    return mass
end

function add_stiff_matrix!(el::Element1d, mass::AbstractMatrix)
end

mass_matrix(el::Element1d) = jacobian(el) * mass_matrix(basis1d(el))
stiff_matrix(el::Element1d) = jacobian(el)^3 * stiff_matrix(basis1d(el))



project(el::Element1d, f::AbstractVector) = project(basis1d(el), f)
project(el::Element1d, f::Function) = project(basis1d(el), qnodes(el))
