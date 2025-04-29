export higher_order_symbolics!,higher_order_state

function higher_order_symbolics!(symbolic_timeseries::Vector{T}, order::Int; delays = collect(0:order-1),warn_overflow=true) where T<:Integer
    length(delays) != order && throw(ArgumentError("The number of provided delays has to be the same as `order`"))
    delays[1] != 0 && throw(ArgumentError("The first delay has to be zero, viz. the actual element."))
    overflow=false
    nr_symbols = maximum(symbolic_timeseries) #assumed symbols go from 1 to nr_symbols

    #create new symbolic timeseries by OVERWRITING the old timeseries vector
    order>1 ? max_delay=maximum(delays) : max_delay=0

    #helper vector to store symbol indices
    #has order nr of elements
    idxs = zeros(Int,order)

    #loop through symb time series
    #i can take order:length(symbolic_timeseries) in default case
    for i in 1+max_delay:length(symbolic_timeseries)
        
        #store indices of symbols used to build higher order symbol
        # has "order" number of elements
        for d in 1:order
            idxs[d] = i-delays[d]
        end
        #use @view here?
        hs,of = higher_order_state(symbolic_timeseries[idxs], nr_symbols)
        overflow = of
        # place higher order state back in the original vector
        symbolic_timeseries[i-max_delay] = hs
    end
    warn_overflow && overflow && @warn "Overflowing integer symbols when calculating higher order states!"
    return nothing
end


function higher_order_state(symbols::Vector{T}, nr_symbols;warn_overflow=false) where T<:Integer
    overflow=false #bool for warnings in case of overflow
    hs::T = 0
    for (i,s) in enumerate(symbols)
        multiplier = nr_symbols^(i-1)
        #warn if overflow happens
        multiplier < 0 && (overflow = true)
        hs += s*multiplier
    end
    warn_overflow && overflow && @warn "Overflowing integer symbols when calculating higher order states!"
    return hs,overflow
end