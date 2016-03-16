
function add_mass_matrix!{T}(bas::GenBasis1d, el::Element1d, mass::AbstractMatrix{T}, λ=one(T))
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

mass_matrix{T}(bas::GenBasis1d, el::Element1d{T}, λ=one(T)) =
    add_mass_matrix!(bas, el, zeros(T, nmodes(bas), nmodes(bas)), λ)



function add_mass_matrix!{T}(el::Element1d, mass::AbstractMatrix{T}, λ::AbstractVector{T})
    M = nmodes(bas)
    Q = nquad(bas)
    ϕ = basis1d(bas)
    wJ = jacweights(el)
    for k = 1:M
        for i = k:M
            m = 0.0
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


function add_stiff_matrix!(el::Element1d, mass::AbstractMatrix)
end


   
