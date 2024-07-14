export randomwalk_step,random_walk_on_stn,network_measures

function randomwalk_step(P::SparseMatrixCSC,src)
	src_row = P[src,:]
	dst = sample(src_row.nzind, Weights(nonzeros(src_row)))
	
	return dst,-log(P[src,dst])
end

"""
Conducts a random walk process on a STN for `N_steps`.
Returns the normalized walk length.
"""
#took out transient, and uses randomwalk_step now
function random_walk_on_stn(P::SparseMatrixCSC, N_steps)
    src = sample(1:size(P)[1]);
    walk_length = 0.0;

    for _ in 1:N_steps
		src, l = randomwalk_step(P,src)
		walk_length = walk_length + l
    end
    
    return walk_length
end

"""
	network_measures(stn::MetaDiGraph,ensemble,N_steps) -> S, Λ
Calculates the Sinai-Kolmogorov Entropy and Lyapunov measure of a STN
by calculating the average walk length and the variance of walk lenghts
over an ensemble of random walks on a STN
"""
function network_measures(P::SparseMatrixCSC, ensemble, N_steps)
    walk_length = Vector{Float64}(undef, ensemble)
    for i in 1:ensemble
        walk_length[i] = random_walk_on_stn(P, N_steps)
    end
   	entropy = mean(walk_length)/N_steps
    lyapunov_measure = var(walk_length,corrected=false)/N_steps
    return entropy, lyapunov_measure
end



"""
	measure_convergence(stn::MetaDiGraph,ensemble,N_max) -> entropy_timeseries,lyapunov_timeseries
Calculates and returns network measures for an ensemble at every step in the random walk up to N_max.
"""
#calc variance without correction
function measure_convergence(stn,ensemble,N_max)
	error("This function is not updated yet!")
	ensemble_walk_lengths = [] #container for walk lengths for every step for every random_walk
	
	for i in 1:ensemble
		source = sample(1:nv(stn));
		walk_length_timeseries = zeros(N_max) #container for individual walk lengths for every step 
		walk_length = 0
		P = get_transition_matrix(stn)
		
		for n in 1:N_max
			source, l = randomwalk_step(stn,source,P) #make one step in the graph
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
