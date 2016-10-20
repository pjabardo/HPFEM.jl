


" Generic basis "
abstract GenBasis{T<:Number}

" Generic 1d basis "
abstract GenBasis1d{T<:Number} <: GenBasis{T}

"""
Structure that stores 1d Basis function information.

In FEM, integrals and derivatives are calculated in standard elements, 
which in the case of 1D domains is the interval (-1,1). 

This type stores the information obtained from the combination
of a set of basis functions and a set of quadrature rules.
"""
type Basis1d{T<:Number,B<:BasisFun1d} <: GenBasis1d{T}
    "Number of modes"
    M::Int

    "Number of quadrature nodes"
    Q::Int

    "Quadrature nodes"
    ξ::Vector{T}

    "Quadrature weights"
    w::Vector{T}

    "Derivative matrix for functions known at quadrature nodes"
    D::Array{T,2}

    "Basis functions at quadrature nodes"
    ϕ::Array{T,2}

    "Derivative of basis functions at quadrature nodes"
    dϕ::Array{T,2}

    "Weight X Basis"
    wϕ::Array{T,2}
    
    "Inverse of mass matrix"
    imass::Cholesky{T}

    "Actual basis function. Used to compute the basis function"
    bas::B
end

    
"""
Construct a Basis from a set of basis functions and quadrature rules.
"""    
function Basis1d{T<:Number, B<:BasisFun1d}(b::B, q::QuadType, ::Type{T}=Float64)
    # Create the basis function
    m = nmodes(b)
    # Obter as informações da quadratura
    Q = q.Q
    ξ = q.z
    w = q.w
    D = q.D
    ϕ = zeros(T, Q, m)
    wϕ = zeros(T, Q, m)
    
    # Preencher as funções de base:
    
    for k = 1:m
        for i=1:Q
            ϕ[i,k] = basis1d(b, ξ[i], k)
            wϕ[i,k] = w[i] * ϕ[i,k]
        end
    end
    
    # Calcular as derivadas:
    dϕ = D * ϕ
    
    # calcular a matrix de massa
    mass = ϕ' * wϕ #std_mass_matrix1d(ϕ, w)
    
    imass = cholfact(Hermitian(mass))
    Basis1d{T,B}(m, Q, ξ, w, D, ϕ, dϕ, wϕ, imass, b)
end


"Return quadrature nodes"
qnodes(b::Basis1d) = b.ξ

"Return number of modes"
nmodes(b::Basis1d) = b.M

"Return number of quadrature nodes"
nquad(b::Basis1d) = b.Q

"Quadrature weights"
qweights(b::Basis1d) = b.w

"Return basis function"
basis(b::Basis1d) = b.bas

"Basis functions at quadrature nodes"
qbasis(b::Basis1d) = b.ϕ

"Derivative matrix at quadrature nodes"
diffmat(b::Basis1d) = b.D

"Derivative of basis functions at quadrature nodes"
dqbasis(b::Basis1d) = b.dϕ

locmap(b::Basis1d) = b.bas.lnum


basis1d(b::GenBasis1d, x, p) = basis1d(basis(b), x, p)

basis1d!(b::GenBasis1d, x::AbstractArray, y::AbstractArray, p) = basis1d!(basis(b), x, y, p)
basis1d(b::GenBasis1d, x::AbstractArray, p) = basis1d!(b, x, similar(x), p)

#call(b::GenBasis1d, x, p) = basis1d(basis(b), x, p)

Basis1d(m::Int, q::Int) = Basis1d(ModalC01d(m), QuadType(q))
Basis1d(m::Int) = Basis1d(m, m+1)

nbndry(b::GenBasis1d) = nbndry(basis(b))
ninterior(b::GenBasis1d) = ninterior(basis(b))

bndidx(b::GenBasis1d) = bndidx(basis(b))
intidx(b::GenBasis1d) = intidx(basis(b))

seq2bi!{T}(b::GenBasis1d, x::AbstractVector{T}, y::AbstractVector{T}) = seq2bi!(basis(b).lnum, x, y)
seq2bi{T}(b::GenBasis1d, x::AbstractVector{T}) = seq2bi(basis(b).lnum, x)
bi2seq!{T}(b::GenBasis1d, x::AbstractVector{T}, y::AbstractVector{T}) = bi2seq!(basis(b).lnum, x, y)
bi2seq{T}(b::GenBasis1d, x::AbstractVector{T}) = bi2seq(basis(b).lnum, x)
bi2seq!{T}(b::GenBasis1d, xb::AbstractVector{T}, xi::AbstractVector{T}, y::AbstractVector{T}) =
    bi2seq!(basis(b).lnum, xb, xi, y)
bi2seq{T}(b::GenBasis1d, xb::AbstractVector{T}, xi::AbstractVector{T}) = bi2seq(basis(b).lnum, xb, xi)
seq2b!{T}(b::GenBasis1d, x::AbstractVector{T}, y::AbstractVector{T}) = seq2b!(basis(b).lnum, x, y)
seq2i!{T}(b::GenBasis1d, x::AbstractVector{T}, y::AbstractVector{T}) = seq2i!(basis(b).lnum, x, y)
seq2b{T}(b::GenBasis1d, x::AbstractVector{T}) = seq2b(basis(b).lnum, x)
seq2i{T}(b::GenBasis1d, x::AbstractVector{T}) = seq2i(basis(b).lnum, x)


function std_mass_matrix1d{T<:Number}(ϕ::AbstractArray{T,2}, w::AbstractVector{T})
    m = size(ϕ,2)
    Q = length(w)
    mass = zeros(T,m,m)
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
    

    return mass
end

""" Calculate the mass matrix


"""
mass_matrix(b::GenBasis1d) =  std_mass_matrix1d(qbasis(b), qweights(b))


"""
Project a function know at quadrature nodes to a basis
"""
function project!(b::Basis1d, f::AbstractVector, u::AbstractVector)

  ϕ = qbasis(b)
  w = qweights(b)
  Q = nquad(b)
  M = nmodes(b)
  iM = b.imass

  for k = 1:M
    F = zero(eltype(f))
    for q = 1:Q
      F += f[q] * ϕ[q,k] * w[q]
    end
    u[k] = F
  end

  A_ldiv_B!(iM, u)

  return u

end

project(b::Basis1d, f::AbstractVector) = project!(b, f, similar(f, nmodes(b)))
"""
Project a function  to a basis

"""
project(b::Basis1d, f::Function) = project(b, f(qnodes(b)))








type SEM1d{T<:Number} <: GenBasis1d
    "Number of quadrature nodes"
    Q::Int

    "Quadrature nodes"
    ξ::Vector{T}

    "Quadrature weights"
    w::Vector{T}

    "Derivative matrix for functions known at quadrature nodes"
    D::Array{T,2}

    "Actual basis function. Used to compute the basis function"
    bas::Lagrange1d{T}

    
end

function SEM1d{T<:Number}(n::Integer, ::Type{T}=Float64)
    q = QuadType(n, Jacobi.GLJ, T)
    SEM1d{T}(q.Q, q.z, q.w, q.D, Lagrange1d(q.z))
end

qnodes(b::SEM1d) = b.ξ
nmodes(b::SEM1d) = b.Q
nquad(b::SEM1d) = b.Q
qweights(b::SEM1d) = b.w
basis(b::SEM1d) = b.bas
qbasis{T<:Number}(b::SEM1d{T}) = eye(T,nmodes(b))
dqbasis(b::SEM1d) = b.D
diffmat(b::SEM1d) = b.D

locmap(b::SEM1d) = b.bas.lnum

project!(b::SEM1d, f::AbstractVector, u::AbstractVector) = copy!(u, f)
project(b::SEM1d, f::AbstractVector) = project!(b, f, similar(f))
project(b::SEM1d, f::Function) = project(b, f(qnodes(b)))



