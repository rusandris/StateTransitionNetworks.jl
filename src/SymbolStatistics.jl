export SymbolStatistics,SymbolStatisticsNoRemap
export add_transition!,add_transitions!
export calculate_transition_matrix,update_dict!,normalize_rows!,calculate_symbol_probabilities

abstract type AbstractSymbolStatistics end 

#--------------------SymbolStatistics------------------

mutable struct SymbolStatistics{T} <: AbstractSymbolStatistics                    
    symbol_space_size::Int64
    Q::SparseMatrixCOO{Float64, Int64}
    x::SparseVector{Float64,Int64}
    symbol_dict::Dict{T,Int64}
    last_transition::Tuple{T,T}
    nr_unique_symbols::Int64 
    nr_transitions::Int64 
end


function SymbolStatistics(T::Type,symbol_space_size::Int64)
    Q = SparseMatrixCOO{Float64, Int64}(Int64[], Int64[], Float64[],symbol_space_size,symbol_space_size)
    x = SparseVector(symbol_space_size,Int64[],Float64[])
    symbol_dict = Dict{T,Int64}()
    return SymbolStatistics{T}(symbol_space_size,Q,x,symbol_dict,(zero(T),zero(T)),0,0)
end

function add_transition!(trs::SymbolStatistics,tr::Tuple{S,S} where S<:Any)

    #src and dst symbols
    src = tr[1]
    dst = tr[2]

    #update internal symbol_dict
    update_dict!(trs,src)
    update_dict!(trs,dst)
    
    #add transition to sparse matrix (setindex! pushes)
    trs.Q[trs.symbol_dict[src] ,trs.symbol_dict[dst] ] = 1.0

    #count occurrence of symbols in x 
    if src != trs.last_transition[2] 
        trs.x[trs.symbol_dict[src]] += 1.0   
        trs.x[trs.symbol_dict[dst]] += 1.0   
    else
        trs.x[trs.symbol_dict[dst]] += 1.0
    end

    #increment nr_transitions
    trs.nr_transitions += 1

    #save last_transition
    trs.last_transition = tr

    return nothing

end

function update_dict!(trs::SymbolStatistics,symbol::S where S<:Any)
    #check if symbol is already used
    if !haskey(trs.symbol_dict,symbol)
        #check if symbol can be added
        if trs.nr_unique_symbols < trs.symbol_space_size
            #push new symbol_dict entry
            push!(trs.symbol_dict, symbol => trs.nr_unique_symbols + 1)
            trs.nr_unique_symbols += 1
        else
            error("Cannot add more symbols! trs.used_symbols ($(trs.nr_unique_symbols)) < trs.symbol_space_size ($(trs.symbol_space_size)) ")
        end
    end
    return nothing
end


function calculate_transition_matrix(trs::SymbolStatistics)
    #transitions are put in order in trs.Q
    nzrows = 1:trs.nr_unique_symbols
    #easy slice
    P = SparseMatrixCSC(trs.Q)
    P = P[nzrows,nzrows]
    calculate_transition_matrix!(P)
    return P
end                     

#--------------------SymbolStatisticsNoRemap------------------
#NoRemap means symbols and their indices in transition matrix is the same

mutable struct SymbolStatisticsNoRemap{I<: Integer} <: AbstractSymbolStatistics                    
    symbol_space_size::I
    Q::SparseMatrixCOO{Float64, I}
    x::SparseVector{Float64,I}
    last_transition::Tuple{I,I}
    nr_unique_symbols::I 
    nr_transitions::I 
end

function SymbolStatisticsNoRemap(I::Type,symbol_space_size::Int64)
    Q = SparseMatrixCOO{Float64, I}(I[], I[], Float64[],symbol_space_size,symbol_space_size)
    x = SparseVector(symbol_space_size,I[],Float64[])
    return SymbolStatisticsNoRemap(symbol_space_size,Q,x,(zero(I),zero(I)),0,0)
end

function add_transition!(trs::SymbolStatisticsNoRemap,tr::Tuple{I,I} ) where I<:Integer

    #src and dst symbols
    src = tr[1]
    dst = tr[2]
    
    #add transition to sparse matrix (setindex! pushes)
    trs.Q[src,dst] = 1.0

    #count occurrence of symbols in x 
    if src != trs.last_transition[2] 
        trs.x[src] += 1.0
        trs.x[dst] += 1.0
    else
        trs.x[dst] += 1.0
    end

    #increment nr_transitions
    trs.nr_transitions += 1

    #save last_transition
    trs.last_transition = tr

    #update nr_unique_symbols
    trs.nr_unique_symbols = length(trs.x.nzind)

    return nothing

end

function calculate_transition_matrix(trs::SymbolStatisticsNoRemap)
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


function add_transitions!(trs::AbstractSymbolStatistics, symbolic_timeseries)  

    n = length(symbolic_timeseries)

    for i in 1:n-1 
        add_transition!(trs,(symbolic_timeseries[i],symbolic_timeseries[i+1]))
    end

end

function calculate_symbol_probabilities(trs::T) where T <: AbstractSymbolStatistics
    x_nz = nonzeros(trs.x) 
    return x_nz ./ sum(x_nz)  
end