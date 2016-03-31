
function add_mass_matrix!{T<:Number}(bas::Basis1d{T}, el::Element1d, mass::AbstractMatrix{T}, λ=one(T))
    M = nmodes(bas)
    Q = nquad(bas)
    ϕ = qbasis(bas)
    wJ = jacweights(el)
    for k = 1:M
        for i = k:M
            m = 0.0
            for q = 1:Q
                m += ϕ[q,k] * ϕ[q, i] * wJ[q] 
            end
            mass[k,i] += m*λ
            if k != i
                mass[i,k] += m*λ
            end
        end
    end
    return mass
end

mass_matrix{T<:Number}(bas::Basis1d{T}, el::Element1d{T}, λ=one(T)) =
    add_mass_matrix!(bas, el, zeros(T, nmodes(bas), nmodes(bas)), λ)

   
stiff_matrix{T<:Number}(bas::Basis1d{T}, el::Element1d{T}, λ=one(T)) =
    add_stiff_matrix!(bas, el, zeros(T, nmodes(bas), nmodes(bas)), λ)


function add_mass_matrix!{T<:Number}(bas::Basis1d{T}, el::Element1d, mass::AbstractMatrix{T}, λ::AbstractVector{T})
    M = nmodes(bas)
    Q = nquad(bas)
    ϕ = qbasis(bas)
    wJ = jacweights(el)
    for k = 1:M
        for i = k:M
            m = zero(T)
            for q = 1:Q
                m += ϕ[q,k] * ϕ[q, i] * wJ[q] * λ[q]
            end
            mass[k,i] += m
            if k != i
                mass[i,k] += m
            end
        end
    end
    return mass
end


function add_stiff_matrix!{T<:Number}(bas::Basis1d{T}, el::Element1d, mat::AbstractMatrix{T}, λ=one(T))
    M = nmodes(bas)
    Q = nquad(bas)
    ϕ = qbasis(bas)
    dϕ = dqbasis(bas)
    D = diffmat(bas)
    wJ = jacweights(el)
    dξdx = deriv_ξ(el)
    for k = 1:M
        for i = k:M
            L = zero(T)
            for q = 1:Q
                L += dϕ[q,i] * dϕ[q,k] *wJ[q] * (dξdx[q]^2)*λ
            end
            mat[i,k] += L
            if k != i
                mat[k,i] += L
            end
        end
    end
    return mat
end

stiff_matrix{T<:Number}(bas::Basis1d{T}, el::Element1d{T}, λ=one(T)) =
    add_stiff_matrix!(bas, el, zeros(T, nmodes(bas), nmodes(bas)), λ)

function add_stiff_matrix!{T<:Number}(bas::Basis1d{T}, el::Element1d, mat::AbstractMatrix{T},
                                      λ::AbstractVector{T})
    M = nmodes(bas)
    Q = nquad(bas)
    ϕ = qbasis(bas)
    dϕ = dqbasis(bas)
    D = diffmat(bas)
    wJ = jacweights(el)
    dξdx = deriv_ξ(el)
    println("DOIS")

    for k = 1:M
        for i = k:M
            L = zero(T)
            for q = 1:Q
                L += dϕ[q,i] * dϕ[q,k] *wJ[q] * (dξdx[q]^2)*λ[q]
            end
            mat[i,k] += L
            if k != i
                mat[k,i] += L
            end
        end
    end

    return mat
end


   
