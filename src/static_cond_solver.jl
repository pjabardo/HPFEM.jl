abstract LinearSolver

abstract BBSolver
abstract StaticCond <: LinearSolver


type CholeskySC{T<:BBSolver} <: StaticCond
    dof::DofMap
    Abb::T
    iAii::Array{Array{Float64, 2},1}
    Abi_iAii::Array{Array{Float64, 2},1}



using Base.LinAlg.BLAS.gemm!
using Base.LinAlg.LAPACK.potrf!
using Base.LinAlg.LAPACK.potrs!
using ArrayViews

function add_local_matrix(solver::CholeskySC, e, Ae)
    dof = solver.dof
    iAii = solver.iAii[e]
    ni = num_bi(dof, e)
    nb = num_be(dof, e)
    for i = 1:ni
        ii = i + nb
        for k = i:ni
            iAii[k, i] = Ae[k+nb, ii]
        end
    end
    potrf!('L', iAii )


    BiC = solver.Abi_iAii[e]
    for k = 1:nb
        for i = 1:ni
            BiC[i,k] = Ae[i+nb, k]
        end
    end

    potrs!('L', iAii, BiC)

    Abb = view(Ae, 1:nb, 1:nb)

    Aib = view(Ae, (nb+1):(nb+ni), 1:nb)

    gemm!('T', 'N', -1.0, BiC, Aib, 1.0, Abb)

    # Now we should assemble the global boundary-boundary matrix.
    # assemble(solver, Abb, e)
end







