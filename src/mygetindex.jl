function getindex_I_sorted_linear(A::AbstractSparseMatrixCSC{Tv,Ti}, I::AbstractVector, J::AbstractVector) where {Tv,Ti}
    #println("WARNING: Experimental getindex_I_sorted_linear")
    require_one_based_indexing(A, I, J)
    nI = length(I)
    nJ = length(J)

    colptrA = getcolptr(A); rowvalA = rowvals(A); nzvalA = nonzeros(A)
    colptrS = Vector{Ti}(undef, nJ+1)
    colptrS[1] = 1

    m = size(A,1)

    # since we linear in complexity we can just use a hashmap
    hashmapI = Dict(I[i] => i for i = nI:-1:1)

    ptrS   = 1
    # build the cache and determine result size
    @inbounds for j = 1:nJ
        col = J[j]
        ptrI::Int = 1 # runs through I
        ptrA::Int = colptrA[col]
        stopA::Int = colptrA[col+1]
        while ptrI <= nI && ptrA < stopA
            rowA = rowvalA[ptrA]
            rowI = I[ptrI]

            if rowI > rowA
                ptrA += 1
            elseif rowI < rowA
                ptrI += 1
            else # found a match
                haskey(hashmapI, rowA) || (hashmapI[rowA] = ptrI)
                ptrS += 1
                ptrI += 1
            end
        end
        colptrS[j+1] = ptrS
    end

    rowvalS = Vector{Ti}(undef, ptrS-1)
    nzvalS  = Vector{Tv}(undef, ptrS-1)

    # fill the values
    ptrS = 1
    @inbounds for j = 1:nJ
        col = J[j]
        ptrA::Int = colptrA[col]
        stopA::Int = colptrA[col+1]
        while ptrA < stopA
            rowA = rowvalA[ptrA]
            ptrI = haskey(hashmapI, rowA) ? hashmapI[rowA] : 0
            if ptrI > 0
                while ptrI <= nI && I[ptrI] == rowA
                    rowvalS[ptrS] = ptrI
                    nzvalS[ptrS] = nzvalA[ptrA]
                    ptrS += 1
                    ptrI += 1
                end
            end
            ptrA += 1
        end
    end
    return SparseMatrixCSC(nI, nJ, colptrS, rowvalS, nzvalS)
end

function getindex_I_sorted_bsearch_I(A::AbstractSparseMatrixCSC{Tv,Ti}, I::AbstractVector, J::AbstractVector) where {Tv,Ti}
    #println("WARNING: Experimental getindex_I_sorted_bsearch_I")
    require_one_based_indexing(A, I, J)
    nI = length(I)
    nJ = length(J)

    colptrA = getcolptr(A); rowvalA = rowvals(A); nzvalA = nonzeros(A)
    colptrS = Vector{Ti}(undef, nJ+1)
    colptrS[1] = 1

    m = size(A, 1)

    # scout the length of cache we need
    ptrS::Int = 0
    @inbounds for j = 1:nJ
        col = J[j]
        ptrS += (colptrA[col+1] - colptrA[col])
    end

    # new definition of cacheI uses the same indexing as rowval of A[:, J] (in compressed format),
    # so indexing into it is easy, while avoiding the complexity depending on m
    # Additionally we use a hashmap to avoid repeatedly searching for the same entry
    cacheI = zeros(Int, ptrS)

    # find relevant rows and put them into the 
    ptrS = 0 # this will count the number of nonzeros that we will have in S
    ptrI::Int = 1 # this is the position in the cache
    ptrC::Int = 1 # points to the current position in the cache
    @inbounds for j = 1:nJ
        col = J[j]
        ptrA = colptrA[col]
        stopA = colptrA[col+1]
        while ptrA < stopA
            rowA = rowvalA[ptrA]
            ptrI = searchsortedfirst(I, rowA, ptrI, nI, Base.Order.Forward)
            @inbounds cacheI[ptrC] = ptrI
            while ptrI <= nI && I[ptrI] == rowA
                ptrS += 1
                ptrI += 1
            end
            (ptrI > nI) && break
            ptrA += 1
            ptrC += 1
        end
        ptrC += stopA - ptrA # adjust cache pointer
        ptrI = 1
    end
    rowvalS = Vector{Ti}(undef, ptrS)
    nzvalS  = Vector{Tv}(undef, ptrS)
    colptrS[nJ+1] = ptrS+1

    # second pass to fill the values (this one doesn't perform bsearch so it is much faster)
    ptrS = 1
    ptrC = 1
    @inbounds for j = 1:nJ
        col = J[j]
        ptrA = colptrA[col]
        stopA = colptrA[col+1]
        while ptrA < stopA 
            rowA = rowvalA[ptrA]
            ptrI = cacheI[ptrC]
            (ptrI > nI) && break
            while ptrI > 0 && I[ptrI] == rowA
                @inbounds rowvalS[ptrS] = ptrI
                @inbounds nzvalS[ptrS] = nzvalA[ptrA]
                ptrS += 1
                ptrI += 1
                (ptrI > nI) && break
            end
            ptrA += 1
            ptrC += 1
        end
        ptrC += stopA - ptrA
        @inbounds colptrS[j+1] = ptrS
    end
    return SparseMatrixCSC(nI, nJ, colptrS, rowvalS, nzvalS)
end