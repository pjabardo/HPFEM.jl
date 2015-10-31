abstract AbstractElement1d
    
type Element1d
  basis::Basis1d
  a::Float64
  b::Float64
  J::Float64
  Element1d(bas::Basis1d, a=-1.0, b=1.0) = new(bas, a, b, 2/(b-a))
end

num_modes(el::Element1d) = el.basis.M
num_quad(el::Element1d) = el.basis.Q
basis1d(el::Element1d) = el.basis
basis(el::Element1d) = basis(basis1d(el))
dbasis(el::Element1d) = dbasis(basis1d(el))
weights(el::Element1d) = weights(basis1d(el))
jacobian(el::Element1d) = el.J

function qnodes(el::Element1d)
  ξ = qnodes(basis1d(el))
  (1 - ξ)*el.a/2 + (1 + ξ)*el.b/2
end


mass_matrix(el::Element1d) = jacobian(el) * mass_matrix(basis1d(el))
stiff_matrix(el::Element1d) = jacobian(el)^3 * stiff_matrix(basis1d(el))



project(el::Element1d, f::AbstractVector) = project(basis1d(el), f)
project(el::Element1d, f::Function) = project(basis1d(el), qnodes(el))
