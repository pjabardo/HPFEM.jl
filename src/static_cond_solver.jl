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
    ub::Vector{T}
    lft::Dict{Int,DirichiletLift}
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

    lmap = locmap(dof)

    for i = 1:nel
        nbe = nbndry(lmap)
        nie = niinterior(lmap)
        Aii[i] = zeros(T, nie, nie)
        M[i] = zeros(T, nie, nbe)
    end
    ub = zeros(T, nbslv)
    lft = Dict{Int,DirichiletLift}()
    CholeskySC(dof, Abb, Aii, M, ub, lft)
    
end

using Base.LinAlg.BLAS.gemm!
using Base.LinAlg.LAPACK.potrf!
using Base.LinAlg.LAPACK.potrs!

function add_local_matrix{Mat<:BBSolver, T<:Number}(solver::CholeskySC{T, Mat}, e::Integer,
                                                    Ae::AbstractMatrix{T}, ib, ii)
    lmap = locmap(solver.dof)
    nb = nbndry(lmap)
    ni = ninterior(lmap)

    if hasdirbc(solver.dof, e)
        lft[e] = DirichiletLift(Ae, idirbc[e])
    end

    Aii = solver.Aii[e]
    for i = 1:ni
        for k = 1:ni
            Aii[k,i] = Ae[ ii[k], ii[i] ]
        end
    end
    potrf!('L', Aii)

    M = solver.M[e]
    for k = 1:nb
        for i = 1:ni
            M[i,k] = Abi[ib[k],ii[i]]
        end
    end

    potrs!('L', Aii, M)
    ib = bndry_idx(lmap)
    ii = interior_idx(lmap)

    Abb = Ae[ib,ib]
    Abi = Ae[ib,ii]
    gemm!('N', 'N', -1.0, Abi, M, 1.0, Abb)

    
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


function add_local_rhs{Mat<:BBSolver, T<:Number}(solver::CholeskySC{Mat,T}, e::Integer,
                                                 Fb::AbstractVector{T}, Fi::AbstractVector{T})
    if hasdirbc(solver.dof, e)
        lift!(solver.lft[e], Fb)
    end
    

end



type LU_SC{Mat<:BBSolver, T <: Number} <: StaticCond
    dof::DofMap
    Abb::Mat
    Aii::Array{T, 3}
    M::Array{T,3}
end






