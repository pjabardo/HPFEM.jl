"Abstract base class for linear system solvers"
abstract LinearSolver


"Abstract class for solvers using static condensation"
abstract StaticCond <: LinearSolver


"Solver using Static condensation on symmetric positive definite problems."
type CholeskySC{Mat<:BBSolver, Dof <: DofMap, T <: Number} <: StaticCond
    dof::Dof
    Abb::Mat
    Aii::Array{T, 3}
    M::Array{T,3}
end


using Base.LinAlg.BLAS.gemm!
using Base.LinAlg.LAPACK.potrf!
using Base.LinAlg.LAPACK.potrs!
using ArrayViews

function add_local_matrix{Mat<:BBSolver, T <: Number}(solver::CholeskySC{Mat,T}, e, Ae)
    dof = solver.dof
    Aii = view(solver.Aii, :, :)
    ni = num_bi(dof, e)
    nb = num_be(dof, e)
    for i = 1:ni
        ii = i + nb
        for k = i:ni
            Aii[k, i] = Ae[k+nb, ii]
        end
    end
    potrf!('L', Aii )


    M = view(solver.M, :, :, e)
    for k = 1:nb
        for i = 1:ni
            M[i,k] = Ae[i+nb, k]
        end
    end

    potrs!('L', Aii, M)
    Abb = zeros(T, nb, nb)

    Aib = view(Ae, (nb+1):(nb+ni), 1:nb)

    gemm!('T', 'N', -1.0, M, Aib, 1.0, Abb)

    # Now we should assemble the global boundary-boundary matrix.
    assemble(solver, Abb, e)
end


type LU_SC{Mat<:BBSolver, T <: Number} <: StaticCond
    dof::DofMap
    Abb::Mat
    Aii::Array{T, 3}
    M::Array{T,3}
end






