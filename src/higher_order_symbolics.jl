export higher_order_symbolics!

function higher_order_symbolics!(symbolic_timeseries::Vector{T}, order::Int; delays = collect(0:order-1)) where T<:Integer
    length(delays) != order && @error "The number of provided delays has to be the same as `order`"
    delays[1] != 0 && @error "The first delay has to be zero, viz. the actual element."

    nr_symbols = maximum(symbolic_timeseries)
    # new symbolic timeseries by OVERWRITING the old timeseries vector
    order>1 ? max_delay=maximum(delays) : max_delay=0
    for i in 1+max_delay:length(symbolic_timeseries)
        # build the higher order vector state
        hs::T = 0
        for d in 1:order
            s = symbolic_timeseries[i-delays[d]]
            hs += s*nr_symbols^(d-1)
        end
        # encode the higher order state
        symbolic_timeseries[i-max_delay] = hs
    end
    return nothing
end


function higher_order_state(states::Vector{T}, nr_symbols::T) where T<:Integer
    hs::T = 0
    for (i,s) in enumerate(states)
        hs += s*nr_symbols^(i-1)
    end
    return hs
end