# DofMap
abstract DofMap

type DofMap <: DofMap
  nb::Int
  nblsv::Int
  nbdir::Int
  nel::Int
  nbe::Int
  nie::Int
  nloc::Int
  bmap::Array{Int,2}
  map::Array{Int,2}

  function DofMap(nel, nb, nbdir, nbe, nie, bmap::Array{Int,2})
    nbslv = nb - nbdir
    nloc = nie + nbe
    mp = zeros(Int, nloc, nel)
    for e = 1:nel
      mp[1,e] = bmap[1,e]
      mp[2,e] = bmap[2,e]
      for i = 1:nie
        mp[2+i,e] = nb + (e-1)*nie + i
      end
    end
    new(nb, nbslv, nbdir, nel, nbe, nie, nloc, bmap, mp)
  end
end


export DofMap

function DofMap(M, nnodes, idir, iper=false)
  nbe = 2
  if M > 2
    nie = M-2
  else
    nie = 0
  end

  nel = nnodes - 1
  bmap = zeros(Int, nbe, nel)

  for e = 1:nel
    bmap[1,e] = e
    bmap[2,e] = e+1
  end
  nb = nnodes
  nd = length(idir)
  if iper
    bmap[2,nel] = 1
    nb = nnodes - 1
    nd = 0
  elseif nd > 0 # There is a Dirichilet BC.
    if nd == 1
      if idir[1] == 1
        ii = [nnodes;  1:(nnodes-1)]
      else
        ii = [1:nnodes]
      end
    else
      ii = [nnodes; 1:(nnodes-1)]
    end

    for e = 1:nel
      bmap[1,e] = ii[e]
      bmap[2,e] = ii[e+1]
    end
  end

  return DofMap(nel, nb, nd, nbe, nie, bmap)
end

num_be(dof::DofMap, e) = dof.nbe
num_ie(dof::DofMap, e) = dof.nbi

function global2local(dof::DofMap, xg::Array{Float64,1})

  xloc = zeros(dof.nloc, dof.nel)
  for e = 1:dof.nel
    xloc[:,e] = xg[dof.map[:,e]]
  end
  return xloc
end

