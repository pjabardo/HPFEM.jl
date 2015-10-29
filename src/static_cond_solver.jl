"Abstract base class for linear system solvers"
abstract LinearSolver


"Abstract class for solvers using static condensation"
abstract StaticCond <: LinearSolver


"Solver using Static condensation on symmetric positive definite problems."
type CholeskySC{Mat<:BBSolver, Dof <: DofMap, T <: Number} <: StaticCond
    dof::Dof
    Abb::Mat
    Aii::Vector{Matrix{T}}
    M::Vector{Matrix{T}}
end


using Base.LinAlg.BLAS.gemm!
using Base.LinAlg.LAPACK.potrf!
using Base.LinAlg.LAPACK.potrs!

function add_local_matrix{Mat<:BBSolver, T <: Number}(solver::CholeskySC{Mat,T},
                                                      e::Integer, Ae::Array{T,2},
                                                      mape, mapi)
    dof = solver.dof
    Aii = solver.Aii[e]
    ni = length(mape)
    nb = length(mapi)
    for i = 1:ni
        for k = i:ni
            Aii[k, i] = Ae[mapi[k], mapi[i]]
        end
    end
    potrf!('L', Aii )


    M = solver.M[e]
    for k = 1:nb
        for i = 1:ni
            M[i,k] = Ae[mapi[i], mape[k]] #Ae[i+nb, k]
        end
    end

    potrs!('L', Aii, M)
    Abb = zeros(T, nb, nb)

    Aib = zeros(T, ni, nb)  #(Ae, (nb+1):(nb+ni), 1:nb)
    for k = 1:nb
        for i = 1:ni
            Aib[i,k] = Ae[mapi[i], mape[k]]
        end
    end
            

    gemm!('T', 'N', -1.0, M, Aib, 1.0, Abb)

    # Now we should assemble the global boundary-boundary matrix.
    assemble(solver, Abb, e, 1:nb)
end


type LU_SC{Mat<:BBSolver, T <: Number} <: StaticCond
    dof::DofMap
    Abb::Mat
    Aii::Array{T, 3}
    M::Array{T,3}
end






