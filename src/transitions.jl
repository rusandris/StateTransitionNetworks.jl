export IntegerTransitions,GeneralTransitions
export add_transition!,add_transitions!
export calculate_transition_matrix,update_dict!,normalize_rows!,calculate_symbol_probabilities

abstract type AbstractTransitions end 

#--------------------IntegerTransitions------------------

mutable struct IntegerTransitions{I<: Integer} <: AbstractTransitions                    
    symbol_space_size::I
    Q::SparseMatrixCOO{Float64, I}
    x::SparseVector{Float64,I}
    nr_transitions::I
end

function IntegerTransitions(I::Type,symbol_space_size::Int64)
    Q = SparseMatrixCOO{Float64, I}(I[], I[], Float64[],symbol_space_size,symbol_space_size)
    x = SparseVector(symbol_space_size,I[],Float64[])
    nr_transitions::I = 0
    return IntegerTransitions(symbol_space_size,Q,x,nr_transitions)
end

function add_transition!(trs::IntegerTransitions,tr::Tuple{I,I} ) where I<:Integer

    #src and dst symbols
    src = tr[1]
    dst = tr[2]
    
    #add transition to sparse matrix (setindex! pushes)
    trs.Q[src,dst] = 1.0

    #count occurrence of symbols in x 
    if trs.nr_transitions == 0
        trs.x[src] += 1.0
        trs.x[dst] += 1.0
    else
        trs.x[dst] += 1.0
    end

    #increment nr_transitions
    trs.nr_transitions += 1

    return nothing

end

function calculate_transition_matrix(trs::IntegerTransitions)
    #get used symbols from occupied indices of the sparse COO
    #trs.Q.is is the row coord, trs.Q.js is the col coord of elements
    used_symbols = trs.x.nzind

    #convert from COO to CSC
    P = SparseMatrixCSC(trs.Q)
    #get row and column of every used symbol
    P = P[used_symbols,used_symbols] 
    #slice empty rows and columns
    calculate_transition_matrix!(P)
    return P
end


function add_transitions!(trs::AbstractTransitions, symbolic_timeseries)  

    n = length(symbolic_timeseries)

    for i in 1:n-1 
        add_transition!(trs,(symbolic_timeseries[i],symbolic_timeseries[i+1]))
    end

end

function calculate_symbol_probabilities(trs::T) where T <: AbstractTransitions
    x_nz = nonzeros(trs.x) 
    return x_nz ./ sum(x_nz)  
end

#--------------------GeneralTransitions------------------

mutable struct GeneralTransitions{S} <: AbstractTransitions                    
    symbol_space_size::Int64
    Q::SparseMatrixCOO{Float64, Int64}
    x::SparseVector{Float64,Int64}
    symbol_dict::Dict{S,Int64}
    nr_used_symbols::Int64 
    nr_transitions::Int64 
end


function GeneralTransitions(S::Type,symbol_space_size::Int64)
    Q = SparseMatrixCOO{Float64, Int64}(Int64[], Int64[], Float64[],symbol_space_size,symbol_space_size)
    x = SparseVector(symbol_space_size,Int64[],Float64[])
    symbol_dict = Dict{S,Int64}()
    return GeneralTransitions{S}(symbol_space_size,Q,x,symbol_dict,0,0)
end

function add_transition!(trs::GeneralTransitions,tr::Tuple{S,S} where S<:Any)
    #src and dst symbols
    src = tr[1]
    dst = tr[2]

    #update internal symbol_dict
    update_dict!(trs,src)
    update_dict!(trs,dst)
    
    #add transition to sparse matrix (setindex! pushes)
    trs.Q[trs.symbol_dict[src] ,trs.symbol_dict[dst] ] = 1.0

    #increment nr_transitions
    trs.nr_transitions += 1

    #count occurrence of symbols in x 
    if trs.nr_transitions == 0
        trs.x[trs.symbol_dict[src]] += 1.0   
        trs.x[trs.symbol_dict[dst]] += 1.0   
    else
        trs.x[trs.symbol_dict[dst]] += 1.0
    end

    return nothing

end

function update_dict!(trs::GeneralTransitions,symbol::S where S<:Any)
    #check if symbol is already used
    if !haskey(trs.symbol_dict,symbol)
        #check if symbol can be added
        if trs.nr_used_symbols < trs.symbol_space_size
            #push new symbol_dict entry
            push!(trs.symbol_dict, symbol => trs.nr_used_symbols + 1)
            trs.nr_used_symbols += 1
        else
            error("Cannot add more symbols! trs.used_symbols ($(trs.nr_used_symbols)) < trs.symbol_space_size ($(trs.symbol_space_size)) ")
        end
    end
    return nothing
end


function calculate_transition_matrix(trs::GeneralTransitions)
    #transitions are put in order in trs.Q
    nzrows = 1:trs.nr_used_symbols
    #easy slice
    P = SparseMatrixCSC(trs.Q)
    P = P[nzrows,nzrows]
    calculate_transition_matrix!(P)
    return P
end

