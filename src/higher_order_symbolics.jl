export higher_order_symbolics!

function higher_order_symbolics!(symbolic_timeseries, order; delays = collect(0:order-1))
    length(delays) != order && @error "The number of provided delays has to be the same as `order`"
    delays[1] != 0 && @error "The first delay has to be zero, viz. the actual element."
    # rename symbols from 1 to nr_symbols
    symbols = sort!(unique(symbolic_timeseries))
	nr_symbols = length(symbols)
    symbol_dictionary = Dict(symbols .=> 1:nr_symbols)
    # higher order state in vector formalism
    states = zeros(order)
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
    return nothing
end
