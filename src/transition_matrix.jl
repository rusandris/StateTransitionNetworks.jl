export calculate_transition_matrix,calculate_transition_matrix!,calculate_transition_matrix_no_remap,is_stochastic,is_strongly_connected,make_strongly_connected,strongly_connected_components

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
function calculate_transition_matrix(symbolic_timeseries::AbstractArray{T};symbol_dictionary=Dict{T,Int64}(),verbose=true) where T <: Integer
	symbolic_timeseries_copy = deepcopy(symbolic_timeseries)
	P,Q,state_probabilities = calculate_transition_matrix!(symbolic_timeseries_copy;symbol_dictionary=symbol_dictionary,verbose=verbose)
	return P,Q,state_probabilities
end

function calculate_transition_matrix!(symbolic_timeseries::AbstractArray{T};symbol_dictionary=Dict{T,Int64}(),verbose=true) where T <: Integer
	L = length(symbolic_timeseries)
    
	if isempty(symbol_dictionary)
		symbols = unique(symbolic_timeseries)
		nr_symbols = length(symbols)
		symbol_dictionary = Dict(symbols .=> 1:nr_symbols)
	end

	#weight and transition probability matrices
	nr_symbols = length(keys(symbol_dictionary))
	Q = spzeros(nr_symbols, nr_symbols)

	#probability distribution of states
	state_probabilities = zeros(nr_symbols) 

	state = symbol_dictionary[symbolic_timeseries[1]]
	symbolic_timeseries[1] = state
	next_state = 0
	#count transitions and map symbols to indices
	for i in 1:(L - 1)
		next_state = symbol_dictionary[symbolic_timeseries[i+1]]
		Q[state,next_state] += 1
		state_probabilities[state] += 1
		symbolic_timeseries[i+1] = next_state
		state = next_state 
	end
	state_probabilities[symbolic_timeseries[end]] += 1.0
	
    #normalize state distribution
	state_probabilities = state_probabilities ./ sum(state_probabilities)

	#normalize Q and fill P by normalizing rows
    Q .= Q./sum(Q)
	P = calculate_transition_matrix(Q; verbose=verbose)
	
	return P,Q,state_probabilities

end

#method that doesn't use a dictionary
function calculate_transition_matrix_no_remap(symbolic_timeseries;verbose=true) 
	L = length(symbolic_timeseries)

	#assuming symbols are integers 1:n
	#assuming indices are the symbols themselves 
	nr_symbols = maximum(symbolic_timeseries)

	#weight and transition probability matrices
	Q = spzeros(nr_symbols, nr_symbols)

	#probability distribution of states
	state_probabilities = zeros(nr_symbols) 

	#count transitions, assuming indices are the symbols themselves 
	for i in 1:(L - 1)
		Q[symbolic_timeseries[i],symbolic_timeseries[i+1]] += 1
		state_probabilities[symbolic_timeseries[i]] += 1.0
	end
	state_probabilities[symbolic_timeseries[end]] += 1.0

	 #normalize state distribution
	 state_probabilities = state_probabilities ./ sum(state_probabilities)

	 #normalize Q and fill P by normalizing rows
	 Q .= Q./sum(Q)
	 P = calculate_transition_matrix(Q; verbose=verbose)

	 return P,Q,state_probabilities

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

function make_strongly_connected(P::SparseMatrixCSC)
	strong_comps = strongly_connected_components(P) #get strongly connected components
	largest_strong_comp = strong_comps[argmax(length.(strong_comps))] #largest
	sort!(largest_strong_comp)
	nc = length(largest_strong_comp) #length of largest -> Pc : nc x nc matrix

	Pc = spzeros((nc,nc))
	rows = rowvals(P)
	vals = nonzeros(P)
	n,_ = size(P)

	#loop through P
	#select rows and columns only from largest_strong_comp
	for j in 1:length(largest_strong_comp)
		n = largest_strong_comp[j] #jth state from largest_strong_comp
		vals_from_col = nzrange(P, n) #all indices from col
		for v in vals_from_col
			i = rows[v] #row in P

			#skip rows that aren't in largest_strong_comp
			if !(i in largest_strong_comp)
				continue
			end

			ic = findfirst(x -> x == i, largest_strong_comp) #row in Pc
			#write val into Pc
			val = vals[v]
			Pc[ic,j] = val
	   end
	end

	#renormalization step
	calculate_transition_matrix!(Pc)
	return Pc

end

is_strongly_connected(P::SparseMatrixCSC) = is_strongly_connected(SimpleDiGraph(P))
strongly_connected_components(P::SparseMatrixCSC) = strongly_connected_components(SimpleDiGraph(P)) 