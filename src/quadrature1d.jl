

type QuadType{T<:Number}
    Q::Int
    z::Vector{T}
    w::Vector{T}
    D::Matrix{T}
end
function QuadType{T<:Number,QT<:QUADRATURE_TYPE}(q::Int, ::Type{QT}=GLJ, ::Type{T}=Float64, a=0, b=0)
    z = qzeros(QT, q, a, b, T)
    w = qweights(QT, z, a, b)
    D = qdiff(QT, z, a, b)
    QuadType(q, z, w, D)
end
