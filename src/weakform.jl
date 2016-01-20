abstract BilinearForm
abstract LinearForm


"""
Represents a weak form of a differential equation.

To provide flexibility, different operators can be used.
This type is used to store operators of the following form:

$$
\sum_{i=1}^N A_i(u,v) = \sum_{k=1}^M b_k(v)
$$

"""
type WeakForm{Dof<:DofMap}
    A::Vector{BilinearForm}
    Ab::Dict{Int,Vector{BilinearForm}}
    F::Vector{LinearForm}
    Fb::Dict{Int,Vector{LinearForm}}
    dlift::Dict{Int,DirichiletLift}
end


