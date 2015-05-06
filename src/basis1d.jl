using Jacobi

abstract BasisFun
abstract BasisFun1d
import Base.length

immutable ModalC0Fun1d
    N::Int
    bi::Bool
    e1::Int
    e2::Int
    ModalC0Fun1d(N) = new(N, true, 1, N)
end
length(b::ModalC0Fun1d) = b.N
function order(b::ModalC0Fun1d, p)
    if p==1
        o = 1
    elseif p==2
        o = 1
    else
        o = p
    end
    o
end
order(b::ModalC0Fun1d) = N

function basis(b::ModalC0Fun1d, p, x)
    if p==1
        y = (one(x) - x) / 2
    elseif p==length(b)
        y = (one(x) + x) / 2
    else
        y = (one(x) - x) * (one(x) + x) * jacobi(x, p-2, one(x), one(x))/4
    end
end

function basis!(b::ModalC0Fun1d, p, x::AbstractArray, y::AbstractArray)
    for i = 1:length(x)
        y[i] = basis(b, p, x[i])
    end
    y
end

basis(b::ModalC0Fun1d, p, x::AbstractArray) = basis!(b, p, x, zeros(x))


type Basis1d{BF,T}
    N::Int
    Q::Int
    z::Array{1,T}
    w::Array{1,T}
    D::Array{2,T}
    B::Array{2,T}

    function Basis1d{BF<:BasisFun, T<:FloatingPoint}(b::BF, quad::Quadrature{T})

        N = length(b)
    end
end




    


