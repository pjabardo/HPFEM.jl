"Abstract base class for linear system solvers"
abstract LinearSolver


"Abstract class for solvers using static condensation"
abstract StaticCond <: LinearSolver


"Solver using Static condensation on symmetric positive definite problems."
type CholeskySC{T <: Number, Mat<:BBSolver, Dof <: DofMap} <: StaticCond
    dof::Dof
    Abb::Mat
    Aii::Vector{Matrix{T}}
    M::Vector{Matrix{T}}
end

bbmatrix(solver::CholeskySC) = solver.Abb
dofmap(solver::CholeskySC) = solver.dof

function CholeskySC{T<:Number, Mat<:BBSolver, Dof <: DofMap}(dof::Dof, ::Type{Mat}, 
                                                             ::Type{T}=Float64)
    nel = num_elems(dof)
    nbslv  = nbslvmodes(dof)
    nb = nbmodes(dof)

    Abb = Mat{T}(nb, nbslv)
    Aii = Vector{Array{T,2}}(nel)
    M = Vector{Array{T,2}}(nel)

    for i = 1:nel
        nbe = nbemodes(dof, i)
        nie = niemodes(dof, i)
        Aii[i] = zeros(T, nie, nie)
        M[i] = zeros(T, nie, nbe)
    end

    CholeskySC(dof, Abb, Aii, M)
    
end

using Base.LinAlg.BLAS.gemm!
using Base.LinAlg.LAPACK.potrf!
using Base.LinAlg.LAPACK.potrs!

function add_local_matrix{Mat<:BBSolver, T<:Number}(solver::CholeskySC{T, Mat}, e::Integer,
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
    
    assemble!(bbmatrix(solver), Abb, bmap(dofmap(solver), e))
end


    
    
function teste_add_local_matrix{Mat<:BBSolver, T <: Number}(solver::CholeskySC{Mat,T},
                                                      e::Integer, Ae::Array{T,2},
                                                      ib::AbstractVector{Int},
                                                      ii::AbstractVector{Int})
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






