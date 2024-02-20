export IntegerTransitions,GeneralTransitions
export add_transition!,add_transitions!
export calculate_transition_matrix,update_dict!,normalize_rows!

#--------------------IntegerTransitions------------------

mutable struct IntegerTransitions{I<: Integer}                     
    symbol_space_size::I
    Q::SparseMatrixCOO{Float64, I}
end

function IntegerTransitions(I::Type,symbol_space_size::Int64)
    Q = SparseMatrixCOO{Float64, I}(I[], I[], Float64[],symbol_space_size,symbol_space_size)
    return IntegerTransitions(symbol_space_size,Q)
end

function add_transition!(trs::IntegerTransitions,tr::Tuple{I,I} ) where I<:Integer

    #src and dst symbols
    src = tr[1]
    dst = tr[2]
    
    #add transition to sparse matrix (setindex! pushes)
    trs.Q[src,dst] = 1.0

    return nothing

end

function calculate_transition_matrix(trs::IntegerTransitions)
    #get used symbols from occupied indices of the sparse COO
    #trs.Q.is is the row coord, trs.Q.js is the col coord of elements
    used_symbols = unique(reduce(vcat,[trs.Q.is,trs.Q.js]))

    #convert from COO to CSC
    P = sparse(trs.Q)
    #get row and column of every used symbol
    P = P[used_symbols,used_symbols] 
    #slice empty rows and columns
    normalize_rows!(P)
    return P
end


function add_transitions!(trs::IntegerTransitions,symbolic_timeseries::Vector{I}) where I <: Integer

    n = length(symbolic_timeseries)

    for i in 1:n-1 
        add_transition!(trs,(symbolic_timeseries[i],symbolic_timeseries[i+1]))
    end

end



#--------------------GeneralTransitions------------------

mutable struct GeneralTransitions{S}                    
    symbol_space_size::Int64
    Q::SparseMatrixCOO{Float64, Int64}
    symbol_dict::Dict{S,Int64}
    nr_used_symbols::Int64 
end


function GeneralTransitions(symbol_space_size::Int64,S::Type)
    Q = SparseMatrixCOO{Float64, Int64}(Int64[], Int64[], Float64[],symbol_space_size,symbol_space_size)
    symbol_dict = Dict{S,Int64}()
    return GeneralTransitions{S}(symbol_space_size,Q,symbol_dict,0)
end

function add_transition!(trs::GeneralTransitions,tr::Tuple{S,S} where S<:Any)
    #src and dst symbols
    src = tr[1]
    dst = tr[2]

    #update internal symbol_dict
    update_dict!(trs,src)
    update_dict!(trs,dst)
    
    #add transition to sparse matrix
    trs.Q[trs.symbol_dict[src] ,trs.symbol_dict[dst] ] += 1.0

    return nothing

end

function update_dict!(trs::GeneralTransitions,symbol::S where S<:Any)
    #check if symbol is already used
    if !haskey(trs.symbol_dict,symbol)
        #check if symbol can be added
        if trs.used_symbols < trs.symbol_space_size
            #push new symbol_dict entry
            push!(trs.symbol_dict, symbol => trs.used_symbols + 1)
            trs.used_symbols += 1
        else
            error("Cannot add more symbols! trs.used_symbols ($(trs.used_symbols)) < trs.symbol_space_size ($(trs.symbol_space_size)) ")
        end
    end
    return nothing
end


function calculate_transition_matrix(trs::GeneralTransitions)
    #transitions are put in order in trs.Q
    nzrows = 1:trs.nr_used_symbols
    #easy slice
    P = sparse(trs.Q)
    P = P[nzrows,nzrows]
    normalize_rows!(P)
    return P
end

function normalize_rows!(S::SparseMatrixCSC;verbose=true)

    stochasticity = true

    St = spzeros(size(S))
    ftranspose!(St,S, x -> x)
    vals = nonzeros(St)
    _,n = size(St)

    #loop over columns
	for j in 1:n
        sumSi = 0.0
        #loop nonzero values from that column
        nzi = nzrange(St,j)
        for i in nzi
            sumSi += vals[i]
        end

        #catch rows (columns) with only zero values
        sumSi == 0.0 && (stochasticity=false)

        #normalize
        for i in nzi
            vals[i] /= sumSi
        end
    end
    ftranspose!(S,St, x->x)
    (stochasticity == false && verbose) && @warn "Transition matrix is not stochastic!"
    nothing
end