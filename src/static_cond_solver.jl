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

function add_local_matrix{Mat<:BBSolver, T<:Number}(solver::CholeskySC{Mat,T}, e::Integer,
                                                    Abb::AbstractMatrix{T},
                                                    Abi::AbstractMatrix{T},
                                                    Aii::AbstractMatrix{T})
    nb = size(Abb,1)
    ni = size(Aii,1)

    iAii = solver.Aii[e]
    for i = 1:ni
        for k = 1:ni
            iAii[k,i] = Aii[k,i]
        end
    end
    potrf!('L', iAii)

    M = solver.M[e]
    for k = 1:nb
        for i = 1:ni
            M[i,k] = Abi[k,i]
        end
    end

    potrs!('L', iAii, M)

    gemm!('T', 'T', -1.0, M, Abi, 1.0, Abb)
    
    assemble(solver, Abb, e, bmap(e))
end


    
    
function add_local_matrix{Mat<:BBSolver, T <: Number}(solver::CholeskySC{Mat,T},
                                                      e::Integer, Ae::Array{T,2},
                                                      ib, ii)
    Abb = Ae[ib,ib]
    Aii = Ae[ii,ii]
    Abi = Ae[ib,ii]

    add_local_matrix(solver, e, Abb, Abi, Aii)
end


type LU_SC{Mat<:BBSolver, T <: Number} <: StaticCond
    dof::DofMap
    Abb::Mat
    Aii::Array{T, 3}
    M::Array{T,3}
end






