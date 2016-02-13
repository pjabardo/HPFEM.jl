abstract GenFunction
abstract GenFunction1d

immutable ConstFun1d{T<:Number} <:GenFunction1d
    val::T
end

call(f::ConstFun1d, x) = f.val


immutable Fun1d{Callable, T<:Number} <: GenFunction1d
    fun::Callable
end

call(f::Fun1d, x) = f(x)


