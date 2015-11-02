"""
# Stores local DOF numbering systems

Each local degree of freedom (DOF) has an associated index. 
Often, there is an advantage to separating the degrees of 
freedom into boundary dof (nonzero support on the boundaries
of the element) and interior degrees of freedom (zero on 
boundaries). Depending on how the basis DOF are initially
defined, two local numbering systems are necessary:

 * Sequential numbering, most natural when computing
 * Boundary/Interior numbering where boundary modes are numbered first and then interior modes.

"""
immutable LocalNumSys
    "Number of boundary modes"
    nb::Int
    "Number of interior modes"
    ni::Int
    "Indicesof boundary modes"
    bndry::Vector{Int}
    "Indices of interior modes"
    interior::Vector{Int}
end
LocalNumSys(b::Vector{Int}, i::Vector{Int}) = LocalNumSys(length(b), length(i), b, i)

"Total number of modes"
nmodes(n::LocalNumSys) = n.nb + n.ni

"Number of boundary modes"
nbndry(n::LocalNumSys) = n.nb

"Number of interior modes"
ninterior(n::LocalNumSys) = n.ni

"Indices of boundary modes"
bndry_idx(n::LocalNumSys) = n.bndry

"Indices of interior modes"
interior_idx(n::LocalNumSys) = n.interior

"Convert vector from sequential numbering to boundary/interior"
function seq2bi!{T}(lnum::LocalNumSys, x::AbstractVector{T}, y::AbstractVector{T})
    
    for i = 1:lnum.nb
        y[i] = x[lnum.bndry[i]]
    end
    for i = 1:lnum.ni
        y[i+lnum.nb] = x[lnum.interior[i]]
    end
    y
end

"Convert vector from sequential numbering to boundary/interior"
seq2bi{T}(lnum::LocalNumSys, x::AbstractVector{T}) = seq2bi!(lnum, x, zeros(T, nmodes(lnum)))

"Convert vector from boundary/interior numbering to sequential numbering"
function bi2seq!{T}(lnum::LocalNumSys, x::AbstractVector{T}, y::AbstractVector{T})
    for i = 1:lnum.nb
        y[lnum.bndry[i]] = x[i]
    end
    for i = 1:lnum.ni
        y[lnum.interior[i]] = x[i+lnum.nb]
    end
    y
end

"Convert vector from boundary/interior numbering to sequential numbering"
bi2seq{T}(lnum::LocalNumSys, x::AbstractVector{T}) = bi2seq!(lnum, x, zeros(T, nmodes(lnum)))


"Convert vector from boundary/interior numbering to sequential numbering"
function bi2seq!{T}(lnum::LocalNumSys, xb::AbstractVector{T}, xi::AbstractVector{T},
                    y::AbstractVector{T})
    for i = 1:lnum.nb
        y[lnum.bndry[i]] = xb[i]
    end
    for i = 1:lnum.ni
        y[lnum.interior[i]] = xi[i]
    end
    y
end

"Convert vector from boundary/interior numbering to sequential numbering"
bi2seq{T}(lnum::LocalNumSys, xb::AbstractVector{T}, xi::AbstractVector{T}) =
    bi2seq!(lnum, xb, xi, zeros(T, nmodes(lnum)))


"Extract boundary modes from modes numbered sequentially"
function seq2b!{T}(lnum::LocalNumSys, x::AbstractVector{T}, y::AbstractVector{T})
    for i = 1:lnum.nb
        y[i] = x[lnum.bndry[i]]
    end
    y
end

"Extract boundary modes from modes numbered sequentially"
seq2b{T}(lnum::LocalNumSys, x::AbstractVector{T}) = seq2b!(lnum, x, zeros(T, nbndry(lnum)))


"Extract interior modes from modes numbered sequentially"
function seq2i!{T}(lnum::LocalNumSys, x::AbstractVector{T}, y::AbstractVector{T})
    for i = 1:lnum.ni
        y[i] = x[lnum.interior[i]]
    end
    y
end

"Extract interior modes from modes numbered sequentially"
seq2i{T}(lnum::LocalNumSys, x::AbstractVector{T}) = seq2i!(lnum, x, zeros(T, ninterior(lnum)))




