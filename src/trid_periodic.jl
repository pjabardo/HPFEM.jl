# Implements solver for periodic tridiagonal systems.

type TridiagonalP{T} <: AbstractMatrix{T}
    dl::Vector{T}
    d::Vector{T}
    du::Vector{T}
    du2::Vector{T}
    x2::Vector{T}
end
