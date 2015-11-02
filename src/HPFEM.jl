module HPFEM

# package code goes here
include("localnum.jl")
include("basis1d.jl")
include("element1d.jl")
include("dofmap.jl")
include("dirichilet_lift.jl")
include("bbmatrix.jl")
include("static_cond_solver.jl")
end # module
