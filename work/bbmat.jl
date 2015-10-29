

# Just a work on the

Ae = [3.0 1.0; 1.0 3.0]
Ae2 = [3.0 1.0; -1.0 3.0]
function test_sym_non_per(Ae, nel, nb, nbslv)

    m = zeros(Int, 2, nel)
    for e in 1:nel
        m[1, e] = e
        m[2, e] = e+1
    end
    
    Ag = HPFEM.BBMatrix{Float64}(nb, nbslv)
    
    # Assemble the matrix:
    for e = 1:nel
        HPFEM.assemble!(Ag, Ae, m[:,e])
    end
    A = copy(Ag.A)
    
    HPFEM.trf!(Ag)
    
    
    f = rand(nbslv)
    
    x0 = A\copy(f)
    
    x1 = HPFEM.trs!(Ag, copy(f))
    
    
    Ag_st = HPFEM.BBSymTri{Float64}(nb, nbslv)
    for e = 1:nel
        HPFEM.assemble!(Ag_st, Ae, m[:,e])
    end
    
    HPFEM.trf!(Ag_st)
    x2 = HPFEM.trs!(Ag_st, copy(f))
    
    
    Ag_tri = HPFEM.BBTri{Float64}(nb, nbslv)
    
    for e = 1:nel
        HPFEM.assemble!(Ag_tri, Ae, m[:,e])
    end
    
    HPFEM.trf!(Ag_tri)
    x3 = HPFEM.trs!(Ag_tri, copy(f))
    
    Ag_tp = HPFEM.BBTriP{Float64}(nb, nbslv, false)
    
    for e = 1:nel
        HPFEM.assemble!(Ag_tp, Ae, m[:,e])
    end
    HPFEM.trf!(Ag_tp)
    x4 = HPFEM.trs!(Ag_tp, copy(f))

    hcat(x0, x1-x0, x2-x0, x3-x0, x4-x0)
end





function test_non_per(Ae, nel, nb, nbslv)

    m = zeros(Int, 2, nel)
    for e in 1:nel
        m[1, e] = e
        m[2, e] = e+1
    end
    
    Ag = HPFEM.BBMatrix{Float64}(nb, nbslv)
    
    # Assemble the matrix:
    for e = 1:nel
        HPFEM.assemble!(Ag, Ae, m[:,e])
    end
    A = copy(Ag.A)
    
    HPFEM.trf!(Ag)
    
    
    f = rand(nbslv)
    
    x0 = A\copy(f)
    
    x1 = HPFEM.trs!(Ag, copy(f))
    
    
    
    Ag_tri = HPFEM.BBTri{Float64}(nb, nbslv)
    
    for e = 1:nel
        HPFEM.assemble!(Ag_tri, Ae, m[:,e])
    end
    
    HPFEM.trf!(Ag_tri)
    x2 = HPFEM.trs!(Ag_tri, copy(f))
    
    Ag_tp = HPFEM.BBTriP{Float64}(nb, nbslv, false)
    
    for e = 1:nel
        HPFEM.assemble!(Ag_tp, Ae, m[:,e])
    end
    HPFEM.trf!(Ag_tp)
    x3 = HPFEM.trs!(Ag_tp, copy(f))

    hcat(x0, x1-x0, x2-x0, x3-x0)
end



function test_per(Ae, nel)

    m = zeros(Int, 2, nel)
    for e in 1:nel
        m[1, e] = e
        m[2, e] = e+1
    end
    m[2,nel] = 1
    nb = nel
    nbslv = nb
    Ag = HPFEM.BBMatrix{Float64}(nb, nbslv)
    
    # Assemble the matrix:
    for e = 1:nel
        HPFEM.assemble!(Ag, Ae, m[:,e])
    end
    A = copy(Ag.A)

    HPFEM.trf!(Ag)
    
    
    f = rand(nbslv)
    
    x0 = A\copy(f)
    
    x1 = HPFEM.trs!(Ag, copy(f))
    

    Ag_tp = HPFEM.BBTriP{Float64}(nb, nbslv, true)
    
    for e = 1:nel
        HPFEM.assemble!(Ag_tp, Ae, m[:,e])
    end
    HPFEM.trf!(Ag_tp)
    x2 = HPFEM.trs!(Ag_tp, copy(f))

    hcat(x0, x1-x0, x2-x0)
end

