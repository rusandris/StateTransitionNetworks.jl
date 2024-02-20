export higher_order_symbolics!

function higher_order_symbolics!(symbolic_timeseries::Vector{T}, order; delays = collect(0:order-1)) where T<:Integer
    length(delays) != order && @error "The number of provided delays has to be the same as `order`"
    delays[1] != 0 && @error "The first delay has to be zero, viz. the actual element."
    # rename symbols from 1 to nr_symbols
    symbols::Vector{T} = sort!(unique(symbolic_timeseries))
	nr_symbols::T = length(symbols)
    symbol_dictionary = Dict(symbols .=> collect(UnitRange{T}(1, nr_symbols)))
    # higher order state in vector formalism
    states::Vector{T} = Vector{T}(undef, order)
    # new symbolic timeseries by OVERWRITING the old timeseries vector
    order>1 ? max_delay=maximum(delays) : max_delay=0
    for i in 1+max_delay:length(symbolic_timeseries)
        # build the higher order vector state
        for d in 1:length(states)
            states[d] = symbol_dictionary[symbolic_timeseries[i-delays[d]]]
        end
        # encode the higher order state
        symbolic_timeseries[i-max_delay] = higher_order_state(states, nr_symbols)
    end
    #return @view symbolic_timeseries[1:length(symbolic_timeseries)-max_delay]
    return nothing
end


function higher_order_state(states::Vector{T}, nr_symbols::T) where T<:Integer
    hs::T = 0
    for (i,s) in enumerate(states)
        hs += s*nr_symbols^(i-1)
    end
    return hs
end