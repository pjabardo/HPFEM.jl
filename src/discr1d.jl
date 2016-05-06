abstract GenFunction
abstract GenFunction1d

abstract H1Space

type H1Space1d{T<:Number, B<:GenBasis1d} <: H1Space
    bas::B{T}
    nnodes::Int
    nel::Int
    nodes::Vector{T}
    elems::Vector{Element1d}
end

function H1Space1d{T<:Number, B<:GenBasis1d}(B{T}, nodes::Vector{T})
    nnodes = length(nodes)
    nel = nnodes - 1
    
end


immutable Constant1d{T<:Number} <:GenFunction1d
    val::T
end

call(f::Constant1d, x) = f.val
getindex(f::Constant1d, x, i) = f.val


immutable Fun1d{Callable, T<:Number} <: GenFunction1d
    fun::Callable
end

call(f::Fun1d, x) = f(x)


