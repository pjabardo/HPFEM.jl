module HPFEM

using Jacobi
using Base.LinAlg
using Polynomials

import Base.length
import Base.call
import Base.size

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
include("plotmatrix.jl")

export isbndryint, ModalC01d, Legendre1d, Lagrange1d
export QuadType
export nbndry, ninterior, bndry_idx, interior_idx
export seq2bi!, seq2bi, bi2seq!, bi2seq, seq2b!, seq2b, seq2i!, seq2i
export Basis1d, qnodes, nmodes, nquad, qweights, basis, qbasis, dqbasis, basis1d, basis1d!

end # module
