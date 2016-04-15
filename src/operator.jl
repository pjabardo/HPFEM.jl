
function add_mass_matrix!{T<:Number}(bas::GenBasis1d, el::Element1d, mass::AbstractMatrix{T}, λ=one(T))
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

function add_mass_matrix!{T<:Number}(bas::SEM1d{T}, el::Element1d, mass::AbstractMatrix{T}, λ=one(T))

    M = nmodes(bas)
    wJ = jacweights(el)

    for i = 1:M
        mass[i,i] += wJ[i]*λ
    end
end


mass_matrix{T<:Number}(bas::GenBasis1d, el::Element1d{T}, λ=one(T)) =
    add_mass_matrix!(bas, el, zeros(T, nmodes(bas), nmodes(bas)), λ)

   

function add_mass_matrix!{T<:Number}(bas::SEM1d{T}, el::Element1d, mass::AbstractMatrix{T}, λ::AbstractVector{T})

    M = nmodes(bas)
    wJ = jacweights(el)

    for i = 1:M
        mass[i,i] += wJ[i]*λ[i]
    end
end

function add_mass_matrix!{T<:Number}(bas::GenBasis1d, el::Element1d, mass::AbstractMatrix{T}, λ::AbstractVector{T})
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


function add_stiff_matrix!{T<:Number}(bas::GenBasis1d, el::Element1d, mat::AbstractMatrix{T}, λ=one(T))
    M = nmodes(bas)
    Q = nquad(bas)
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

stiff_matrix{T<:Number}(bas::GenBasis1d, el::Element1d{T}, λ=one(T)) =
    add_stiff_matrix!(bas, el, zeros(T, nmodes(bas), nmodes(bas)), λ)

function add_stiff_matrix!{T<:Number}(bas::GenBasis1d, el::Element1d, mat::AbstractMatrix{T},
                                      λ::AbstractVector{T})
    M = nmodes(bas)
    Q = nquad(bas)
    dϕ = dqbasis(bas)
    D = diffmat(bas)
    wJ = jacweights(el)
    dξdx = deriv_ξ(el)
    

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


# RHS

function add_rhs!{T<:Number}(bas::GenBasis1d, el::Element1d, f::AbstractVector{T},
                            Fe::AbstractVector{T})

    wJ = jacweights(el)
    ϕ = qbasis(bas)
    M = nmodes(bas)
    Q = nquad(bas)

    for k = 1:M
        F = zero(T)
        for q = 1:Q
            F += f[q] * ϕ[q,k] * wJ[q]
        end
        Fe[k] += F
    end
    return Fe 
end


function add_rhs!{T<:Number}(bas::SEM1d{T}, el::Element1d, f::AbstractVector{T},
                            Fe::AbstractVector{T})

    wJ = jacweights(el)
    M = nmodes(bas)

    for k = 1:M
        Fe[k] = f[k] * wJ[k]
    end
    return Fe 
end

