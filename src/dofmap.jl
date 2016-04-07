# DofMap
abstract DofMap



type DofMap1d <: DofMap
    nb::Int
    nbslv::Int
    nbdir::Int
    nel::Int
    lmap::LocalNumSys1d
    bmap::Array{Int,2}
    map::Array{Int,2}
    idir::Dict{Int,Vector{Int}}
    function DofMap1d(nel, nb, nbdir, lmap::LocalNumSys1d,
                      bmap::Array{Int,2}, idir::Dict{Int,Vector{Int}})
        nbslv = nb - nbdir
        nloc = nmodes(lmap)
        mp = zeros(Int, nloc, nel)
        for e = 1:nel
            mp[1,e] = bmap[1,e]
            mp[2,e] = bmap[2,e]
            for i = 1:ninterior(lmap)
                mp[2+i,e] = nb + (e-1)*ninterior(lmap) + i
            end
        end
        new(nb, nbslv, nbdir, nel, lmap, bmap, mp, idir)
    end
end
nbmodes(dof::DofMap) = dof.nb
nbslvmodes(dof::DofMap) = dof.nbslv
num_elems(dof::DofMap) = dof.nel
ninodes(dof::DofMap1d) = dof.nel * dof.nie

locmap(dof::DofMap1d, e=1) = dof.lmap

hasdirbc(dof::DofMap1d, e) = haskey(dof.idir, e)
idirbc(dof::DofMap1d, e) = dof.idir[e]

export DofMap1d

function DofMap1d(lmap, nnodes, idir, iper=false)
    nbe = nbndry(lmap)
    nie = ninterior(lmap)
    
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
                ii = [nnodes;  1:(nnodes-1);]
            else
                ii = [1:nnodes;]
            end
        else
            ii = [nnodes; 1:(nnodes-1);]
        end
        
        for e = 1:nel
            bmap[1,e] = ii[e]
            bmap[2,e] = ii[e+1]
        end
        
    end
    idir = Dict{Int,Vector{Int}}()
    nbslv = nb - nd

    ib = bndry_idx(lmap)
    for e = 1:nel
        if bmap[1,e] > nbslv && bmap[2,e] > nbslv
            idir[e] = [1,2] 
        elseif bmap[1,e] > nbslv
            idir[e] =  [1] 
        elseif bmap[2,e] > nbslv
            idir[e] = [2] 
        end
    end
    return DofMap1d(nel, nb, nd, lmap, bmap, idir)
end

num_be(dof::DofMap, e) = dof.nbe
num_ie(dof::DofMap, e) = dof.nie

bmap(dof, e) = sub(dof.bmap, :, e)

function global2local(dof::DofMap, xg::Array{Float64,1})

    xloc = zeros(dof.nloc, dof.nel)
    for e = 1:dof.nel
        xloc[:,e] = xg[dof.map[:,e]]
    end
    return xloc
end


