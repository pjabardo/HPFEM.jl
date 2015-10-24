"Abstract type that deals with global boundary-boundary matrices"
abstract BBSolver

"Global boundary-boundary Symmetric Positive-Definite matrix"
abstract BBSym <: BBSolver

"Global boundary-boundary symmetric tridiagonal matrices"
type BBTriSym{T<:Number} <: BBSym
    "Number of boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Diagonal elements of the matrix"
    D::Array{T,1}
    "Sub-diagonal elements"
    E::Array{T,1}
    function BBTriSym(nb, nbslv)
        D = zeros(T,nbslv)
        E = zeros(T, nbslv-1)
        new(nb, nbslv, D, E)
    end
end

using Base.LinAlg.LAPACK.pttrf!
using Base.LinAlg.LAPACK.pttrs!

"LU (or Choleksy) decomposition of boundary-boundary system"
trf!(Ag::BBTriSym) = pttrf!(Ag.D, Ag.E)

"Solving linear system using LU (or Cholesky) decomposition"
trs!(Ag::BBTriSym, x) = pttrs!(Ag.D, Ag.E, x)

"""
Assembles the global boundary-boundary matrix.
"""
function assemble!{T}(Ag::BBTriSym{T}, Ae, m)
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


type BBTriP{T<:Number} <: BBSolver
    "Number of independent boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Is the system periodic"
    periodic::Bool
    "Lower diagonal"
    Dl::Vector{T}
    "Main Diagonal"
    D::Vector{T}
    "Upper diagonal"
    Du::Vector{T}
    "Auxiliary memory"
    x2::Vector{T}
    an1::T
    b1::T
    bn1::T
    cn::T
    cn1::T
    Du2::Vector{T}
    ipiv::Vector{T}
    function BBTriP(nb, nbslv, periodic=false)
        if !periodic
            D = zeros(T,nbslv)
            Dl = zeros(T, nbslv-1)
            Du = zeros(T, nbslv-1)
            new(nb, nbslv, D, Dl, Du)
        else
            nbslv = nb-1
            D = zeros(T,nb-1)
            Dl = zeros(T,nb-2)
            Du = zeros(T,nb-2)
            x2 = zeros(T,nb-1)
            zt = zero(T)
            new(nb, nbslv, periodic, Dl, D, Du, x2, zt, zt, zt, zt, zt)
        end
    end
            
end
                    
function assemble!{T}(Ag::BBTriP{T}, Ae, m)
    per = Ag.periodic

    m1 = m[1]
    m2 = m[2]
    
    nb = Ag.nb
    if per
        nbslv = Ag.nb-1
    else
        nbslv = Ag.nbslv
    end    

    ig = m1
    kg = m2

    D = Ag.D
    Dl = Ag.Dl
    Du = Ag.Du

    if !per || (m2 < nb && m2 != 1) # Not periodic or not the last element (if periodic)

        np1 = Ag.nbslv + 1
        
        
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
    elseif m2 == nb # It is periodic AND it is the next to last
        Ag.an1   += Ae[2,2]
        Ag.cn    += Ae[1,2]
        Ag.bn1   += Ae[2,1]
        D[nbslv] += Ae[1,1]
    elseif m2 == 1
        Ag.b1    += Ae[2,1]
        Ag.cn1   += Ae[1,2]
        D[1]     += Ae[2,2]
        Ag.an1   += Ae[1,1]
    end
    
end
                    

function trf!(Ag::BBTriP)
    Dl, D, Du, Ag.Du2, Ag.ipiv = gttrf!(Ag.Dl, Ag.D, Ag.Du)
    return
end
function trs!(Ag::BBTriP, x)
    if !Ag.periodic
        gttrs!('N', Ag.Dl, Ag.D, Ag.Du, Ag.Du2, Ag.ipiv, x)
    else
        n = Ag.nb
        n1 = n-1
        qn1 = view(x, 1:n1, :)
        nrhs = size(x,2)
        gttrs!('N', Ag.Dl, Ag.D, Ag.Du, Ag.Du2, Ag.ipiv, qn1)
        x2 = Ag.x2
        for i = 1:nrhs
            fill!(x2, 0)
            x2[1] = -Ag.b1
            x2[n1] = -Ag.cn
            gttrs!('N', Ag.Dl, Ag.D, Ag.Du, Ag.Du2, Ag.ipiv, x2)
            x[n,i] = (x[n,i] - Ag.cn1*qn1[1,i] - Ag.bn1*qn1[n1,i]) /
            (Ag.an1 + Ag.cn1*x2[1] + Ag.bn1*x2[n1])
            for k = 1:n1
                x[k,i] += x[n,i]*x2[k]
            end
        end
        x
    end
end

