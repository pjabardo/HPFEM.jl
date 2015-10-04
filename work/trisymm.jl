function rand_matrix(n, symm=true, diag=4.0)
    A = randn(n,n)

    if symm
        A = 0.5*(A+A')
    end
    for i = 1:n
        A[i,i] += diag
    end

    A
end

simple_matrix(n) = [2.0 1.0; 1.0 2.0]

using Base.LinAlg.LAPACK.pttrf!
using Base.LinAlg.LAPACK.pttrs!


function exemplo_tri1(n)
    A = 2.0 * eye(n)
    for i = 1:(n-1)
        A[i,i+1] = 1.0
        A[i+1,i] = 1.0
    end

    D = 2.0 * ones(n)
    E = ones(n-1)

    x1 = randn(n)
    x2 = copy(x1)
    y1 = A\x1
    pttrf!(D, E)
    y2 = pttrs!(D, E, x2)

    maxabs(y1-y2)
end

function assemble!{T<:Number}(Ag::Array{T,2}, Ae, m)
    nloc = size(Ae,1)
    for i = 1:nloc
        ig = m[i]
        for k = 1:nloc
            kg = m[k]
            Ag[kg,ig] += Ae[k,i]
        end
    end
end

function exemplo_tri2(n=5)
    ne = n-1

    Ag = zeros(n,n)
    Atri = HPFEM.BBTriSym{Float64}(n, n)
    for e = 1:ne
        Ae = rand_matrix(2)
        m = [e,e+1]
        HPFEM.assemble!(Atri, Ae, m)
        assemble!(Ag, Ae, m)
    end

    x1 = randn(n)
    x2 = copy(x1)

    y1 = Ag\x1
    HPFEM.trf!(Atri)
    y2 = HPFEM.trs!(Atri, x2)

    maxabs(y2-y1)
end

function exemplo_tri3(n=5)

    ne = n-1

    Ag = zeros(n,n)
    Atri = HPFEM.BBTri{Float64}(n, n)
    for e = 1:ne
        Ae = rand_matrix(2, false)
        m = [e,e+1]
        HPFEM.assemble!(Atri, Ae, m)
        assemble!(Ag, Ae, m)
    end

    x1 = randn(n)
    x2 = copy(x1)

    y1 = Ag\x1
    HPFEM.trf!(Atri)
    y2 = HPFEM.trs!(Atri, x2)

    maxabs(x2 - y1)
    hcat(y1, y2, y2-y1)

end
