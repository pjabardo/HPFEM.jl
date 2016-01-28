module HPFEM

using Jacobi
using Base.LinAlg

# package code goes here
include("localnum.jl")
include("basis_fun1d.jl")
include("quadrature1d.jl")
include("basis1d.jl")
include("element1d.jl")
include("dofmap.jl")
include("dirichilet_lift.jl")
include("bbmatrix.jl")
include("static_cond_solver.jl")
end # module
