"Abstract type that deals with global boundary-boundary matrices"
abstract BBSolver
using Base.LinAlg.LAPACK.BlasInt

type BBMatrix{T<:Number} <: BBSolver
    nb::Int
    nbslv::Int
    A::Array{T,2}
    ipiv::Array{BlasInt,1}
    BBMatrix(nb, nbslv) = new(nb, nbslv, zeros(T, nbslv, nbslv))
end

function assemble!{T<:Number}(Ag::BBMatrix{T}, Ae, m)

    nm = length(m)
    n1 = Ag.nbslv+1
    for i = 1:nm
        ig = m[i]
        if ig < n1
            for k = 1:nm
                kg = m[k]
                if kg < n1
                    Ag.A[kg,ig] += Ae[k,i]
                end
            end
        end
    end
end

using Base.LinAlg.LAPACK.getrf!
using Base.LinAlg.LAPACK.getrs!
function trf!(Ag::BBMatrix)
    A, ipiv, info = getrf!(Ag.A)
    Ag.ipiv = ipiv
end
trs!(Ag::BBMatrix, x) = getrs!('N', Ag.A, Ag.ipiv, x)


"Global boundary-boundary Symmetric Positive-Definite matrix"
abstract BBSym <: BBSolver

"Global boundary-boundary symmetric tridiagonal matrices"
type BBSymTri{T<:Number} <: BBSym
    "Number of boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Diagonal elements of the matrix"
    D::Array{T,1}
    "Sub-diagonal elements"
    E::Array{T,1}

    "Tridiagonal Matrix"
    tri::SymTridiagonal{T}
    "Factorization of symmetric tridiagonal matrix"
    fact:: LDLt{T,SymTridiagonal{T}}
    
    function BBSymTri(nb, nbslv)
        D = zeros(T,nbslv)
        E = zeros(T, max(nbslv-1,0))
        tri = SymTridiagonal(D, E)
        new(nb, nbslv, D, E, tri)
    end
end


"LU (or Choleksy) decomposition of boundary-boundary system"
function trf!(Ag::BBSymTri)
    Ag.fact = ldltfact!(Ag.tri)
end

"Solving linear system using LU (or Cholesky) decomposition"
trs!(Ag::BBSymTri, x) = A_ldiv_B!(Ag.fact, x)

"""
Assembles the global boundary-boundary matrix.
"""
function assemble!{T}(Ag::BBSymTri{T}, Ae, m)
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




"Global boundary-boundary symmetric tridiagonal matrices"
type BBSymTri2{T<:Number} <: BBSym
    "Number of boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Diagonal elements of the matrix"
    D::Array{T,1}
    "Sub-diagonal elements"
    E::Array{T,1}
    function BBSymTri2(nb, nbslv)
        D = zeros(T,nbslv)
        E = zeros(T, max(nbslv-1,0))
        new(nb, nbslv, D, E)
    end
end

using Base.LinAlg.LAPACK.pttrf!
using Base.LinAlg.LAPACK.pttrs!

"LU (or Choleksy) decomposition of boundary-boundary system"
trf!(Ag::BBSymTri2) = pttrf!(Ag.D, Ag.E)

"Solving linear system using LU (or Cholesky) decomposition"
trs!(Ag::BBSymTri2, x) = pttrs!(Ag.D, Ag.E, x)

"""
Assembles the global boundary-boundary matrix.
"""
function assemble!{T}(Ag::BBSymTri2{T}, Ae, m)
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
    "Tridiagonal object"
    tri::Tridiagonal{T}
    fact::Base.LinAlg.LU{T,Tridiagonal{T}}
    """
    Inner constructor. Du2 and ipiv are left undefined,
    first the matrix must be assembled.
    """
    function BBTri(nb, nbslv)
        D = zeros(T,nbslv)
        Dl = zeros(T, max(nbslv-1, 0))
        Du = zeros(T, max(nbslv-1, 0))
        tri = Tridiagonal(Dl, D, Du)
        new(nb, nbslv, D, Dl, Du, tri)
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

function trf!(Ag::BBTri)
    Ag.fact = lufact!(Ag.tri)
    return
end
trs!(Ag::BBTri, x) = A_ldiv_B!(Ag.fact, x)


"Global boundary-boundary tridiagonal matrices"
type BBTri2{T<:Number} <: BBSolver
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
    function BBTri2(nb, nbslv)
        D = zeros(T,nbslv)
        Dl = zeros(T, max(nbslv-1, 0))
        Du = zeros(T, max(nbslv-1, 0))
        new(nb, nbslv, D, Dl, Du)
    end
end


function assemble!{T}(Ag::BBTri2{T}, Ae, m)
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

function trf!(Ag::BBTri2)
    Dl, D, Du, Ag.Du2, Ag.ipiv = gttrf!(Ag.Dl, Ag.D, Ag.Du)
    return
end
trs!(Ag::BBTri2, x) = gttrs!('N', Ag.Dl, Ag.D, Ag.Du, Ag.Du2, Ag.ipiv, x)


type BBTriP{T<:Number} <: BBSolver
    "Number of independent boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
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
    tri::Tridiagonal{T}
    fact::LinAlg.LU{T,Tridiagonal{T}}
    
    function BBTriP(nb, nbslv)
        nbslv = nb-1
        D = zeros(T,nbslv)
        Dl = zeros(T,nbslv-1)
        Du = zeros(T,nbslv-1)
        x2 = zeros(T,nbslv)
        tri = Tridiagonal(Dl, D, Du)
        new(nb, nbslv, Dl, D, Du, x2, zero(T), zero(T), zero(T), zero(T), zero(T), tri)
        
    end
            
end

function assemble!{T}(Ag::BBTriP{T}, Ae, m)

    m1 = m[1]
    m2 = m[2]
    
    nb = Ag.nb
    nbslv = Ag.nb-1

    ig = m1
    kg = m2

    D = Ag.D
    Dl = Ag.Dl
    Du = Ag.Du

    if m2 < nb && m2 != 1 #Not the last element 
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
    elseif m2 == nb # It is the next to last
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
    Ag.fact = lufact!(Ag.tri)
    return
end
function trs!(Ag::BBTriP, x)
    n = Ag.nb
    n1 = n-1
    qn1 = sub(x, 1:n1, :)
    nrhs = size(x,2)
    A_ldiv_B!(Ag.fact, qn1)
    x2 = Ag.x2
    for i = 1:nrhs
        fill!(x2, 0)
        x2[1] = -Ag.b1
        x2[n1] = -Ag.cn
        A_ldiv_B!(Ag.fact, x2)
        x[n,i] = (x[n,i] - Ag.cn1*qn1[1,i] - Ag.bn1*qn1[n1,i]) /
            (Ag.an1 + Ag.cn1*x2[1] + Ag.bn1*x2[n1])
        for k = 1:n1
            x[k,i] += x[n,i]*x2[k]
        end
    end
    x

end



type BBSymTriP{T<:Number} <: BBSym
    "Number of independent boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Main Diagonal"
    D::Vector{T}
    "Upper diagonal"
    Du::Vector{T}
    "Auxiliary memory"
    x2::Vector{T}
    an1::T
    cn::T
    cn1::T
    "Tridiagonal Matrix"
    tri::SymTridiagonal{T}
    "Factorization of symmetric tridiagonal matrix"
    fact:: LDLt{T,SymTridiagonal{T}}
    
    function BBSymTriP(nb, nbslv)
        nbslv = nb-1
        D = zeros(T,nb-1)
        Du = zeros(T,nb-2)
        x2 = zeros(T,nb-1)
        tri = SymTridiagonal(D, Du)
        new(nb, nbslv, D, Du, x2, zero(T), zero(T), zero(T), tri)
    end
end



function assemble!{T<:Number}(Ag::BBSymTriP{T}, Ae, m)

    m1 = m[1]
    m2 = m[2]
    
    nb = Ag.nb
    nbslv = Ag.nb-1

    ig = m1
    kg = m2

    D = Ag.D
    Du = Ag.Du
    
    if m2 < nb && m2 != 1 # Not the last element (if periodic)

        np1 = Ag.nbslv + 1
        
        
        if ig < np1
            D[ig] += Ae[1,1]
            if kg < np1
                D[kg] += Ae[2,2]
                if ig < kg
                    Du[ig] += Ae[1,2]
                elseif kg < ig
                    Du[kg] = Ae[2,1]
                end
            end
        elseif kg < np1
            D[kg] += Ae[2,2]
        end
    elseif m2 == nb # It is  the next to last
        Ag.an1   += Ae[2,2]
        Ag.cn    += Ae[1,2]
        D[nbslv] += Ae[1,1]
    elseif m2 == 1
        Ag.cn1   += Ae[1,2]
        D[1]     += Ae[2,2]
        Ag.an1   += Ae[1,1]
    end
    
end




function trf!(Ag::BBSymTriP)
    Ag.fact = ldltfact!(Ag.tri)
    return
end
function trs!(Ag::BBSymTriP, x)
    n = Ag.nb
    n1 = n-1
    qn1 = sub(x, 1:n1, :)
    nrhs = size(x,2)
    A_ldiv_B!(Ag.fact, qn1)
    x2 = Ag.x2
    for i = 1:nrhs
        fill!(x2, 0)
        x2[1] = -Ag.cn1
        x2[n1] = -Ag.cn
        A_ldiv_B!(Ag.fact, x2)
        x[n,i] = (x[n,i] - Ag.cn1*qn1[1,i] - Ag.cn*qn1[n1,i]) /
            (Ag.an1 + Ag.cn1*x2[1] + Ag.cn*x2[n1])
        for k = 1:n1
            x[k,i] += x[n,i]*x2[k]
        end
    end
    x
end
