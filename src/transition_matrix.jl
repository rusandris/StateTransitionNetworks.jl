export calculate_transition_matrix,is_stochastic

"""
	calculate_transition_matrix(time_discrete_ts::TimeSeries,grid_size::Integer;grid_edges=[],returnQ=false) -> P
Calculates the transition matrix `P` from a time-discrete time series `time_discrete_ts` by dividing the state space into cells (`grid_size`x`grid_size`).
If `returnQ` is set to `true`, the weight matrix `Q` is also returned. Grid edges can be specified explicitly with `grid_edges`, otherwise they're inferred from the data.

"""
function calculate_transition_matrix(time_discrete_ts::TimeSeries,grid_size::Integer;grid_edges=[],returnQ=false,return_state_distribution=false)
	error("Not implemented yet!")
	symbolic_timeseries,vertex_positions = timeseries_to_grid(time_discrete_ts,grid_size;grid_edges=grid_edges)
	calculate_transition_matrix(symbolic_timeseries;symbol_dictionary=vertex_positions,returnQ=returnQ,return_state_distribution=return_state_distribution)
end

"""
	calculate_transition_matrix(symbolic_timeseries::TimeSeries,vertex_positions::AbstractDict; returnQ=false) -> P
Calculates the transition matrix `P` from `symbolic_timeseries`. The symbol dictionary (how symbols map to indices of `P`) can be specified with a Dict `symbol_disctionary`. By default they're indexed in the order of appearance.
If `returnQ` is set to `true`, the weight matrix `Q` is also returned. If `return_state_distribution` is set to `true`, the estimated probability distribution over the states is also returned.

"""
function calculate_transition_matrix(symbolic_timeseries::AbstractVector;map_symbols=true,return_everything=false)
	
    
	if map_symbols
		symbols = unique(symbolic_timeseries)
		nr_symbols = length(symbols)

		#weight and transition probability matrices
		Q = spzeros(nr_symbols, nr_symbols)
		P = spzeros(nr_symbols, nr_symbols)

		#probability distribution of states
		state_probabilities = zeros(nr_symbols) 
	
	
		symbol_dictionary = Dict(symbols .=> 1:nr_symbols)
		
		#count transitions and map symbols to indices
		for i in eachindex(symbolic_timeseries[1:end-1])
			state = symbolic_timeseries[i]
			next_state = symbolic_timeseries[i+1]
		    Q[symbol_dictionary[state],symbol_dictionary[next_state]] += 1
		    state_probabilities[symbol_dictionary[state]] += 1
		end
		
    	end_state = symbolic_timeseries[end]
    	state_probabilities[symbol_dictionary[end_state]] += 1
		
	else
		#assuming symbols are integers 1:n
		#assuming indices are the symbols themselves 
		nr_symbols = maximum(symbolic_timeseries)

		#weight and transition probability matrices
		Q = spzeros(nr_symbols, nr_symbols)
		P = spzeros(nr_symbols, nr_symbols)
		
		#probability distribution of states
		state_probabilities = zeros(nr_symbols) 
	
		#count transitions, assuming indices are the symbols themselves 
		for i in eachindex(symbolic_timeseries[1:end-1])
			state = symbolic_timeseries[i]
			next_state = symbolic_timeseries[i+1]
		    Q[state,next_state] += 1
		    state_probabilities[state] += 1
		end
		end_state = symbolic_timeseries[end]
    	state_probabilities[end_state] += 1
		
	end
	
    #normalize state distribution
	state_probabilities = state_probabilities ./ sum(state_probabilities)

	#normalize Q and fill P by normalizing rows
    Q .= Q./sum(Q)
	P = calculate_transition_matrix(Q)

	if return_everything
		if map_symbols
			return P,Q,state_probabilities,symbol_dictionary
		else
			return P,Q,state_probabilities
		end
	else
		return P
	end
end

function calculate_transition_matrix(S::SparseMatrixCSC;verbose=true)
	S_returned = deepcopy(S)
	calculate_transition_matrix!(S_returned,verbose=verbose)
	return S_returned
end

function calculate_transition_matrix!(S::SparseMatrixCSC;verbose=true)

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

function is_stochastic(S::SparseMatrixCSC)

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

        #catch rows (columns) that don't sum up to 1
        if !isapprox(sumSi,1.0) 
            return false
        end
    end
    return true
end