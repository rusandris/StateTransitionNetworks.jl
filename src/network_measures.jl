function prob_matrix(graph)
	nr_vertices = nv(graph)
	P = spzeros(Float64,(nr_vertices,nr_vertices))
	for edge in collect(edges(graph))
		i = edge.src
		j = edge.dst
		P[i,j] = graph[i,j][1]
	end
	return P
end

function weight_matrix(graph)
	nr_vertices = nv(graph)
	Q = spzeros(Float64,(nr_vertices,nr_vertices))
	for edge in collect(edges(graph))
		i = edge.src
		j = edge.dst
		Q[i,j] = graph[i,j][2]
	end
	return Q

end

function randomwalk_step(graph,source,prob_matrix)
	neigh = outneighbors(graph, source)
	neigh_weights = prob_matrix[source,:]
	destination = sample(neigh, Weights(neigh_weights.nzval))
	w = neigh_weights[destination]
	
	return destination,-log(w)
end

"""
Conducts a random walk process on a weighted graph for `N_steps`.
Returns the normalized walk length.
"""
#took out transient, and uses randomwalk_step now
function random_walk_on_weighted_graph(graph, N_steps)
    source = sample(1:nv(graph));
    walk_length = 0.0;
    P = prob_matrix(graph)

    for n in 1:N_steps
		source, l = randomwalk_step(graph,source,P)
		walk_length = walk_length + l
    end
    
    return walk_length
end

"""
	walk_statistics(ensemble,graph::SimpleWeightedDiGraph,N_steps,transient) -> S, Λ
Calculates the Sinai-Kolmogorov Entropy and Lyapunov measure of a STN
by calculating the average walk length and the variance of walk lenghts
over an ensemble of random walks on weighted graph
"""
function walk_statistics(graph, ensemble, N_steps)
    walk_length = Vector{Float64}(undef, ensemble)
    for i in 1:ensemble
        walk_length[i] = random_walk_on_weighted_graph(graph, N_steps)
    end
   	entropy = mean(walk_length)/N_steps
    lyapunov_measure = var(walk_length,corrected=false)/N_steps
    return entropy, lyapunov_measure
end

"""
	measure_convergence(graph::SimpleWeightedDiGraph,ensemble,N_max) -> entropy_timeseries,lyapunov_timeseries
Calculates and returns network measures for an ensemble at every step in the random walk up to N_max.
"""
#calc variance without correction
function measure_convergence(graph,ensemble,N_max)
	ensemble_walk_lengths = [] #container for walk lengths for every step for every random_walk
	
	for i in 1:ensemble
		source = sample(1:nv(graph));
		walk_length_timeseries = zeros(N_max) #container for individual walk lengths for every step 
		walk_length = 0
		P = prob_matrix(graph)
		
		for n in 1:N_max
			source, l = randomwalk_step(graph,source,P) #make one step in the graph
			walk_length += l
			walk_length_timeseries[n] = walk_length  
		end
		push!(ensemble_walk_lengths,walk_length_timeseries) #save individual walk length timeseries
	end
	ensemble_walk_lengths = hcat(ensemble_walk_lengths...) 
	
	#calculate measures for every step
	steps = 1:N_max
	entropy_timeseries = mean(ensemble_walk_lengths,dims=2) ./ steps
	lyapunov_timeseries = var(ensemble_walk_lengths,dims=2,corrected=false) ./steps 
	
	return entropy_timeseries,lyapunov_timeseries
end

"""
	sinai_kolmogorov_entropy(graph_Q; graph_P=nothing) -> S
Calculates analytically the Sinai-Kolmogorov entropy on a STN. The main input `graph_Q` is a SimpleWeightedDiGraph
object containg the occurence probability of each transition (Q[i,j]). The additional `graph_P` argument is another 
graph containg the normalized transition probabilities (P[i,j]). If this is not provided the function calculates automatically these transition
probabilities.
"""
function sinai_kolmogorov_entropy(graph)
    Q = weight_matrix(graph)
    P = prob_matrix(graph)

    entropy = -sum(Q[Q .!=0] .* log.(P[P .!=0]))
    return entropy
end

function sinai_kolmogorov_entropy(Q,P)
    entropy = -sum(Q[Q .!=0] .* log.(P[P .!=0]))
end

