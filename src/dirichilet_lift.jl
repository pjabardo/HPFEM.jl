

type DirichiletLift
  M::Int
  A::Array{Float64,2}
  isx::Array{Int,1}

  function DirichiletLift(A0, idir)
    M = size(A,1)
    isx = trues(M)
    for i in idir
      isx[i] = false
    end
    A = copy(A0)
    new(M, A, isx)
  end
end

using Base.LinAlg.BLAS

function lift!(lft::DirichiletLift, A::Array{Float64}, b::Vector{Float64}, xd::Vector{Float64})

  isx = lft.isx

  for i = 1:lft.M
    if isx[i]
      xd[i] = 0.0
    end
  end
  gemv!('N', -1.0, lft.A, xd, 1.0, b)
end



