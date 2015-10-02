"Abstract type that deals with global boundary-boundary matrices"
abstract BBSolver

"Global boundary-boundary Symmetric Positive-Definite matrix"
abstract BBSymm <: BBSolver

"Global boundary-boundary symmetric tridiagonal matrices"
type BBTriSymm{T<:Number} <: BBSymm
    "Number of boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Diagonal elements of the matrix"
    D::Array{T,1}
    "Sub-diagonal elements"
    E::Array{T,1}
    function BBTriSymm(nb, nbslv)
        D = zeros(T,nbslv)
        E = zeros(T, nbslv-1)
        new(nb, nbslv, D, E)
    end
end

using Base.LinAlg.LAPACK.pttrf!
using Base.LinAlg.LAPACK.pttrs!

"LU (or Choleksy) decomposition of boundary-boundary system"
trf!(Ag::BBTriSymm) = pttrf!(Ag.D, Ag.E)

"Solving linear system using LU (or Cholesky) decomposition"
trs!(Ag::BBTriSymm, x) = pttrs!(Ag.D, Ag.E, x)

"""
Assembles the global boundary-boundary matrix.
"""
function assemble!{T}(Ag::BBTriSymm{T}, Ae, m)
    np1 = Ag.nbslv + 1
    D = Ag.D
    E = Ag.E

    ig = m[1]
    kg = m[2]

    if ig < np1
        D[ig] += Ae[1,1]
        if kg < np1
            D[kg] += Ae[2,2]
            E[min(ig,kg)] += Ae[1,2]
        end
    elseif kg < np1
        D[kg] += Ae[2,2]
    end

end

using Base.LinAlg.LAPACK.BlasInt

"Global boundary-boundary tridiagonal matrices"
type BBTri{T<:Number} <: BBSolver
    "Number of boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Diagonal elements of the matrix"
    D::Array{T,1}
    "Lower sub-diagonal elements"
    Dl::Array{T,1}
    "Upper sub-diagonal elements"
    Du::Array{T,1}
    "Extra storage needed by Lapack"
    Du2::Array{T,1}
    "Pivoting array"
    ipiv::Array{BlasInt,1}

    """
    Inner constructor. Du2 and ipiv are left undefined,
    first the matrix must be assembled.
    """
    function BBTri(nb, nbslv)
        D = zeros(T,nbslv)
        Dl = zeros(T, nbslv-1)
        Du = zeros(T, nbslv-1)
        new(nb, nbslv, D, Dl, Du)
    end
end


function assemble!{T}(Ag::BBTri{T}, Ae, m)
    np1 = Ag.nbslv + 1
    D = Ag.D
    Dl = Ag.Dl
    Du = Ag.Du


    ig = m[1]
    kg = m[2]

    if ig < np1
        D[ig] += Ae[1,1]
        if kg < np1
            D[kg] += Ae[2,2]
            if ig < kg
                Du[ig] += Ae[1,2]
                Dl[ig] += Ae[2,1]
            elseif kg < ig
                Du[kg] = Ae[2,1]
                Dl[kg] = Ae[1,2]
            end
        end
    elseif kg < np1
        D[kg] += Ae[2,2]
    end

end

using Base.LinAlg.LAPACK.gttrf!
using Base.LinAlg.LAPACK.gttrs!

function trf!(Ag::BBTri)
    Dl, D, Du, Ag.Du2, Ag.ipiv = gttrf!(Ag.Dl, Ag.D, Ag.Du)
    return
end
trs!(Ag::BBTri, x) = gttrs!('N', Ag.Dl, Ag.D, Ag.Du, Ag.Du2, Ag.ipiv, x)

