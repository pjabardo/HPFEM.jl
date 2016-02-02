
"""
Implements Dirichilet lift

In FEM, Neuman boundary conditions are naturally satisfied
(hey, they are also known as natural BCs!) but Dirichilet BCs
are another problem. The simplest solution is to use basis functions
that satisfy the boundary conditions. Now, by design, your solution
satisfies the BCs. This usually is a difficult thing to do. Except
when the solution is zero at Dirichilet BCs. This homogeneous case
is easy with basis that split into boundary modes and interior modes:

Interior modes satisy the homogeneous requirements! Just drop the relevant
boundary modes!

Now we still have the issue of non-homogeneous BCs. That is easy as well,
just split the solution into two parts: homogeneous and dirichilet:


u(x) = u^h(x) + u^d(x)


The weak form will take the following form:


a(u, v) = a(u^h + u^d, v) = a(u^h, v) + a(u^d, v) = f(x)


Since u^d is known (Dirichiulet BCs),


a(u^h, v) = f(x) - a(u^d, v)


The `DirichiletLift` type and associated methods implement the correction of
the right hand side of the equation at a local element.

This implementation is very simple and probably not the most efficient. Just
compute the matrix as if it were a common internal element. Store this matrix
and the indicies that correspond to the Dirichilet modes.

When you want to lift the solution, just set the correct Dirichilet modes
and the `lift!` function will set the other modes to zero and subtract
the contribution.

 The DirichiletLift constructor accepts these arguments

 * `A0` local elemental matrix
 * `idir` local index of dirichilet modes

"""
type DirichiletLift{T <: Number}
    "Order of local matrix"
    M::Int

    "Full operator matrix (copy)"
    A::Array{T,2}

    "Index of unknown variables"
    isx::Array{Int,1}

    """
    """
    function DirichiletLift(A0::AbstractArray{T,2}, idir)
        M = size(A0,1)
        isx = trues(M)
        for i in idir
            isx[i] = false
        end
        A = copy(A0)
        new(M, A, isx)
    end
end

using Base.LinAlg.BLAS

"""
Actually lifts the solution

 * `lft` A `DirichiletLift` object that was previously created.
 * `b` A vector withThe right hand side corresponding to the element
 * `xd` 

"""
function lift!{T<: Number}(lft::DirichiletLift, b::AbstractVector{T},
                           xd::AbstractVector{T})

  isx = lft.isx

  for i = 1:lft.M
    if isx[i]
      xd[i] = zero(T)
    end
  end
  gemv!('N', -1.0, lft.A, xd, 1.0, b)
end



