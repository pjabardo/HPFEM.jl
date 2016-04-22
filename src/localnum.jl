abstract LocalNumSys

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
immutable LocalNumSys1d <: LocalNumSys
    "Number of boundary modes"
    nb::Int
    "Number of interior modes"
    ni::Int
    "Indicesof boundary modes"
    bndry::Vector{Int}
    "Indices of interior modes"
    intrr::Vector{Int}
end
LocalNumSys1d(b::Vector{Int}, i::Vector{Int}) = LocalNumSys1d(length(b), length(i), b, i)

isbndryint(n::LocalNumSys1d) = n.ni > 0

"Total number of modes"
nmodes(n::LocalNumSys1d) = n.nb + n.ni

"Number of boundary modes"
nbndry(n::LocalNumSys1d) = n.nb

"Number of interior modes"
nintrr(n::LocalNumSys1d) = n.ni

"Indices of boundary modes"
bndidx(n::LocalNumSys1d) = n.bndry

"Indices of interior modes"
intidx(n::LocalNumSys1d) = n.intrr

"Convert vector from sequential numbering to boundary/interior"
function seq2bi!{T}(lnum::LocalNumSys1d, x::AbstractVector{T}, y::AbstractVector{T})
    
    for i = 1:lnum.nb
        y[i] = x[lnum.bndry[i]]
    end
    for i = 1:lnum.ni
        y[i+lnum.nb] = x[lnum.intrr[i]]
    end
    y
end

"Convert vector from sequential numbering to boundary/interior"
seq2bi{T}(lnum::LocalNumSys1d, x::AbstractVector{T}) = seq2bi!(lnum, x, zeros(T, nmodes(lnum)))

"Convert vector from boundary/interior numbering to sequential numbering"
function bi2seq!{T}(lnum::LocalNumSys1d, x::AbstractVector{T}, y::AbstractVector{T})
    for i = 1:lnum.nb
        y[lnum.bndry[i]] = x[i]
    end
    for i = 1:lnum.ni
        y[lnum.intrr[i]] = x[i+lnum.nb]
    end
    y
end

"Convert vector from boundary/interior numbering to sequential numbering"
bi2seq{T}(lnum::LocalNumSys1d, x::AbstractVector{T}) = bi2seq!(lnum, x, zeros(T, nmodes(lnum)))


"Convert vector from boundary/interior numbering to sequential numbering"
function bi2seq!{T}(lnum::LocalNumSys1d, xb::AbstractVector{T}, xi::AbstractVector{T},
                    y::AbstractVector{T})
    for i = 1:lnum.nb
        y[lnum.bndry[i]] = xb[i]
    end
    for i = 1:lnum.ni
        y[lnum.intrr[i]] = xi[i]
    end
    y
end

"Convert vector from boundary/interior numbering to sequential numbering"
bi2seq{T}(lnum::LocalNumSys1d, xb::AbstractVector{T}, xi::AbstractVector{T}) =
    bi2seq!(lnum, xb, xi, zeros(T, nmodes(lnum)))


"Extract boundary modes from modes numbered sequentially"
function seq2b!{T}(lnum::LocalNumSys1d, x::AbstractVector{T}, y::AbstractVector{T})
    for i = 1:lnum.nb
        y[i] = x[lnum.bndry[i]]
    end
    y
end

"Extract boundary modes from modes numbered sequentially"
seq2b{T}(lnum::LocalNumSys1d, x::AbstractVector{T}) = seq2b!(lnum, x, zeros(T, nbndry(lnum)))


"Extract interior modes from modes numbered sequentially"
function seq2i!{T}(lnum::LocalNumSys1d, x::AbstractVector{T}, y::AbstractVector{T})
    for i = 1:lnum.ni
        y[i] = x[lnum.intrr[i]]
    end
    y
end

"Extract interior modes from modes numbered sequentially"
seq2i{T}(lnum::LocalNumSys1d, x::AbstractVector{T}) = seq2i!(lnum, x, zeros(T, nintrr(lnum)))




