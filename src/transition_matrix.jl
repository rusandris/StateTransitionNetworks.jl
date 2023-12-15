export calculate_transition_matrix
TimeSeries = Union{AbstractStateSpaceSet,AbstractArray}

"""
	calculate_transition_matrix(time_discrete_ts::TimeSeries,grid_size::Integer;grid_edges=[],returnQ=false) -> P
Calculates the transition matrix `P` from a time-discrete time series `time_discrete_ts` by dividing the state space into cells (`grid_size`x`grid_size`).
If `returnQ` is set to `true`, the weight matrix `Q` is also returned. Grid edges can be specified explicitly with `grid_edges`, otherwise they're inferred from the data.

"""
function calculate_transition_matrix(time_discrete_ts::TimeSeries,grid_size::Integer;grid_edges=[],returnQ=false,return_state_distribution=false)
	symbolic_timeseries,vertex_positions = timeseries_to_grid(time_discrete_ts,grid_size;grid_edges=grid_edges)
	calculate_transition_matrix(symbolic_timeseries;symbol_dictionary=vertex_positions,returnQ=returnQ,return_state_distribution=return_state_distribution)
end

"""
	calculate_transition_matrix(symbolic_timeseries::TimeSeries,vertex_positions::AbstractDict; returnQ=false) -> P
Calculates the transition matrix `P` from `symbolic_timeseries`. The symbol dictionary (how symbols map to indices of `P`) can be specified with a Dict `symbol_disctionary`. By default they're indexed in the order of appearance.
If `returnQ` is set to `true`, the weight matrix `Q` is also returned. If `return_state_distribution` is set to `true`, the estimated probability distribution over the states is also returned.

"""
function calculate_transition_matrix(symbolic_timeseries::TimeSeries;symbol_dictionary=Dict(),returnQ=false,return_state_distribution=false)
	
	symbols = unique(symbolic_timeseries)
	nr_symbols = length(symbols)
	
	if length(symbol_dictionary) == 0
		symbol_dictionary = Dict(symbols .=> 1:nr_symbols)
	end
	
	#weight and transition probability matrices
	Q = spzeros(nr_symbols, nr_symbols)
    P = spzeros(nr_symbols, nr_symbols)
    
    #probability distribution of states
    state_probabilities = zeros(nr_symbols)
    
    #count transitions
    for i in eachindex(symbolic_timeseries[1:end-1])
    	state = symbolic_timeseries[i]
    	next_state = symbolic_timeseries[i+1]
        Q[symbol_dictionary[state],symbol_dictionary[next_state]] += 1
        state_probabilities[symbol_dictionary[state]] += 1
    end
    
    #update end (final) state grid
    end_state = symbolic_timeseries[end]
    state_probabilities[symbol_dictionary[end_state]] += 1
    
    #normalize state distribution
	state_probabilities = state_probabilities ./ sum(state_probabilities)

	#normalize Q and fill P by normalizing rows
    Q = Q./sum(Q)
	P = calculate_transition_matrix(Q)

	if returnQ
		if return_state_distribution
			return P,Q,state_probabilities
		else
			return P,Q
		end
	else
		return P
	end
end

"""
	calculate_transition_matrix(Q::AbstractMatrix;verbose=false) -> P
Calculates the transition matrix `P` from  the weight matrix `Q`. Warns if `P` is not stochastic if `verbose` is `true`.

"""
function calculate_transition_matrix(Q::AbstractMatrix;verbose=false)
	P = spzeros(size(Q))
	for i in 1:size(Q)[1]
		sumQi = sum(Q[i,:])
		P[i,:] = Q[i,:]./sumQi
    end
    !all(p -> isfinite(p),P) && verbose && @warn "The matrix is not stochastic!";
    return P
end
