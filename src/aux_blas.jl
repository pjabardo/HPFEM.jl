
import Base.LinAlg.BLAS.gemv!
import Base.LinAlg.BLAS.gemv
import Base.LinAlg.BLAS.gemm!
import Base.LinAlg.BLAS.gemm


function gemv!{T<:Number}(trans::Char, alpha, A::AbstractMatrix{T},
                          x::AbstractVector{T}, beta, y::AbstractVector{T})

    trans = uppercase(trans)
    o = one(T)
    z = zero(T)
    info = 0
    M = size(A,1)
    N = size(A,2)
    
    if !(trans=='N') && !(trans=='T') && !(trans=='C')
        info = 1
    elseif M < 0
        info = 2
    elseif N < 0
        info = 3
    end

    noconj = false
    if trans=='T'
        noconj = true
    end
    
    
    if M==0 || N==0 || (alpha==0 && beta==1)
        return y
    end

    if trans=='N'
        lenx = N
        leny = M
    else
        lenx = M
        leny = N
    end

    kx = 1
    ky = 1
    
    if beta != 1
        if beta == 0
            for i in 1:leny
                y[i] = zero(T)
            end
        else
            for i = 1:leny
                y[i] = beta*y[i]
            end
        end
    end
    if alpha==0
        return y
    end
    
    if trans=='N'

        for j = 1:N
            if x[j] != z
                temp = alpha * x[j]
                for i = 1:M
                    y[i] +=  temp*A[i,j]
                end
            end
        end

    else

        for j = 1:N
            temp = zero(T)
            if noconj
                for i = 1:M
                    temp += A[i,j]*x[i]
                end
            else
                for i = 1:M
                    temp += conj(A[i,j]) * x[i]
                end
            end
            y[j] += alpha*temp
        end
    end

    return y
    
end

function gemv{T<:Number}(trans::Char, alpha, A::AbstractMatrix{T},x::AbstractVector{T})
    gemv!(trans, one(T), A, x, zero(T), similar(x, T, size(A, (trans=='N'?1:2))))
end






function gemm!{T<:Number}(transA::Char, transB::Char, alpha, A::AbstractMatrix{T}, 
                          B::AbstractMatrix{T}, beta, C::AbstractMatrix{T})
    transA = uppercase(transA)
    transB = uppercase(transB)

    M = size(A, transA=='N'?1:2)  #size(A,1)
    K = size(A, transA=='N'?2:1) #size(A,2)
    Kb = size(B, transB=='N'?1:2)
    N = size(B, transB=='N'?2:1) #size(B,2)
    
    notA = transA == 'N'
    notB = transB == 'N'
    conjA = transA == 'C'
    conjB = transB == 'C'
    
    
    info = 0
    if !notA && !conjA && !(transA=='T')
        info = 1
        throw(ArgumentError("transA should be one of 'N', 'T' or 'C'"))
    elseif !notB && !conjB && !(transB=='T')
        info = 2
        throw(ArgumentError("transB should be one of 'N', 'T' or 'C'"))
    elseif M != size(C,1) || N != size(C,2)
        throw(DimensionMismatch("A has size ($nrowA,$ncolA), B has size ($nrowB,$ncolB), C has size $(size(C))"))
    elseif Kb != K
        throw(DimensionMismatch("op(A) has dimensions ($M,$K) and op(B) has dimensions ($Kb,$N) which are incompatible in matrix multiplication"))
    end
    
        
    if M==0 || N==0 || ((alpha==0 || K==0) && beta==1)
        return C
    end

    if alpha==0
        if beta==0
            for j = 1:N
                for i = 1:M
                    C[i,j] = zero(T)
                end
            end
        else
            for j = 1:N
                for i = 1:M
                    C[i,j] = beta*C[i,j]
                end
            end
        end
        return C
    end

    if notB
        if notA
            for j = 1:N
                if beta==0
                    for i = 1:M
                        C[i,j] = zero(T)
                    end
                elseif beta != 1
                    for i = 1:M
                        C[i,j] = beta*C[i,j]
                    end
                end
                for l = 1:K
                    if B[l,j] != zero(T)
                        temp = alpha * B[l,j]
                        for i = 1:M
                            C[i,j] += temp*A[i,l]
                        end
                    end
                end
            end
        elseif conjA
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += conj(A[l,i]) * B[l,j]
                    end
                    if beta==0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        else
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += A[l,i] * B[l,j]
                    end
                    if beta==0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        end
    elseif notA
        if conjB
            for j = 1:N
                if beta==0
                    for i = 1:M
                        C[i,j] = zero(T)
                    end
                elseif beta!=1
                    for i = 1:M
                        C[i,j] = beta*C[i,j]
                    end
                end
                for l = 1:K
                    if B[j,l] != zero(T)
                        temp = alpha * conj(B[j,l])
                        for i = 1:M
                            C[i,j] += temp*A[i,l]
                        end
                    end
                end
            end
        else
            for j = 1:N
                if beta==0
                    for i = 1:M
                        C[i,j] = zero(T)
                    end
                elseif beta != 1
                    for i = 1:M
                        C[i,j] = beta*C[i,j]
                    end
                end
                for l = 1:K
                    if B[j,l] != zero(T)
                        temp = alpha*B[j,l]
                        for i = 1:M
                            C[i,j] += temp*A[i,l]
                        end
                    end
                end
            end
        end
    elseif conjA
        if conjB
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        #println("$j, $i, $l --- A: $(size(A)) --- B: $(size(B))")
                        temp += conj(A[l,i]) * conj(B[j,l])
                    end
                    if beta == 0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        else
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += conj(A[l,i])*B[j,l]
                    end
                    if beta == 0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        end
    else
        if conjB
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += A[l,i]*conj(B[j,l])
                    end
                    if beta==0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        else
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += A[l,i]*B[j,l]
                    end
                    if beta==0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        end
    end
    return C


end

function gemm{T<:Number}(transA::Char, transB::Char, alpha, A::AbstractMatrix{T}, 
                         B::AbstractMatrix{T})
    gemm!(transA, transB, alpha, A, B, zero(T),
          similar(A, T, size(A, transA=='N'?1:2),
                  size(B, transB=='N'?2:1)))
end

