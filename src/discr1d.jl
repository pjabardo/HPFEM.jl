abstract GenFunction{T<:Number}
abstract GenFunction1d{T<:Number}

abstract Discr{T<:Number}

type Discr1d{T<:Number, B<:GenBasis1d, Dof<:DofMap} <: Discr
    bas::B{T}
    dof::Dof
    nnodes::Int
    nel::Int
    nodes::Vector{T}
    elems::Vector{Element1d}
end

function Discr1d{T<:Number, B<:GenBasis1d}(B{T}, msh::Mesh1d{)
    nnodes = length(nodes)
    nel = nnodes - 1
    
end


immutable Constant1d{T<:Number} <:GenFunction1d
    val::T
end

call(f::Constant1d, x) = f.val
getindex(f::Constant1d, x, i) = f.val


immutable Fun1d{T<:Number, } <: GenFunction1d
    fun::Callable
end

call(f::Fun1d, x) = f(x)


