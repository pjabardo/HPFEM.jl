using PyPlot


function plot_matrix(M, posneg=true,eps=1e-5)

    nrow = size(M,1)
    ncol = size(M,2)

    ntot = nrow * ncol

    A = reshape(M, ntot)

    irow = zeros(Int,ntot)
    icol = zeros(Int,ntot)

    er = maxabs(A)*eps
    
    cnt = 1
    for j = 1:ncol
        for i = 1:nrow
            irow[cnt] = i
            icol[cnt] = j
            cnt = cnt + 1
        end
    end

    plot([0, ncol+1, ncol+1, 0, 0], [0, 0, nrow+1, nrow+1, 0], color="black", linewidth=3)
    for i = 1:ncol
        axvline(i, color="lightgray")
    end
    for i = 1:nrow
        axhline(i, color="lightgray")
    end
    
      
    
    if posneg
        ipos = A .> er
        ineg = A .< -er
        
        plot(icol[ipos], nrow + 1 - irow[ipos], "ro")
        plot(icol[ineg], nrow + 1 - irow[ineg], "bo")
    else
        idx = abs(A) .> er
        plot(icol[idx], nrow + 1 - irow[idx], "o", color="black")
    end
    
    return
end

    
