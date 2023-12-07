"""
	sparse_log(A::AbstractMatrix) -> log(A)
Calculates the piecewise logarithm of a sparse matrix without explicitely converting it to dense matrix.
"""
function sparse_log(A::AbstractMatrix)
	P = copy(A)
	vals = nonzeros(P)
	m, n = size(P)
	for j = 1:n
		for i in nzrange(P, j)
			vals[i] = log(vals[i])
		end
	end
	return P
end

function randomwalk_step(stn,source,P)
	neigh = outneighbors(stn, source)
	neigh_weights = P[source,:]
	destination = sample(neigh, Weights(neigh_weights.nzval))
	w = neigh_weights[destination]
	
	return destination,-log(w)
end

"""
Conducts a random walk process on a STN for `N_steps`.
Returns the normalized walk length.
"""
#took out transient, and uses randomwalk_step now
function random_walk_on_stn(stn, N_steps)
    source = sample(1:nv(stn));
    walk_length = 0.0;
    P = get_transition_matrix(stn)

    for n in 1:N_steps
		source, l = randomwalk_step(stn,source,P)
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
function network_measures(stn, ensemble, N_steps)
    walk_length = Vector{Float64}(undef, ensemble)
    for i in 1:ensemble
        walk_length[i] = random_walk_on_stn(stn, N_steps)
    end
   	entropy = mean(walk_length)/N_steps
    lyapunov_measure = var(walk_length,corrected=false)/N_steps
    return entropy, lyapunov_measure
end

"""
	network_measures(P::AbstractMatrix) -> S, Λ
Calculates the Sinai-Kolmogorov Entropy and Lyapunov measure of a STN
by using the analytical definitions of both quantities
"""
function network_measures(P::AbstractMatrix;x=nothing)
   	entropy, ret_code_entr = sinai_kolmogorov_entropy(P;x=x)
    lyapunov, variance, covariance, ret_code_lyap = lyapunov_measure(P;x=x)
	return entropy, lyapunov
end

"""
	measure_convergence(stn::MetaDiGraph,ensemble,N_max) -> entropy_timeseries,lyapunov_timeseries
Calculates and returns network measures for an ensemble at every step in the random walk up to N_max.
"""
#calc variance without correction
function measure_convergence(stn,ensemble,N_max)
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

"""
	sinai_kolmogorov_entropy(stn::MetaDiGraph) -> S
Calculates analytically the Sinai-Kolmogorov entropy on a STN by extracting the Q weight matrix and P transition probability matrix.
"""
function sinai_kolmogorov_entropy(stn)
    Q = get_weight_matrix(stn)
    P = get_transition_matrix(stn)

    entropy = -sum(Q .* -sparse_log(P))
    return entropy
end

"""
	stationary_distribution(P) -> x
Calculates the stationary probability based on the probability matrix P. Each element of the resulting vector is the probability of finding the system in node i.
"""
function stationary_distribution(P::AbstractMatrix)
	x = nullspace(Matrix(P' - I))
	x = x ./sum(x)
	return x
end

"""
	sinai_kolmogorov_entropy(Q,P) -> S
Calculates analytically the Sinai-Kolmogorov entropy given the Q weight matrix and P transition probability matrix of the STN. 

"""
function sinai_kolmogorov_entropy(Q,P)
    entropy = -sum(Q .* -sparse_log(P))
	return entropy
end

"""
	sinai_kolmogorov_entropy(P) -> S
Calculates analytically the Sinai-Kolmogorov entropy given the P transition probability matrix of the STN. 

"""
function sinai_kolmogorov_entropy(P::AbstractMatrix;x=nothing)
		
	if isnothing(x)
		x = stationary_distribution(P)
	end
	
	v = ones(length(x))

	L = -sparse_log(P)
	L = P.*L
	entropy = (x'*L*v)[1]
	return real(entropy), :Success
end


"""
	lyapunov_measure(P) -> Λ
Calculates analytically the Lyapunov measure given the the P transition probability matrix of the STN. 

"""
function lyapunov_measure(P::AbstractMatrix;x=nothing)
	
	if isnothing(x)
		x = stationary_distribution(P)
	end
	x = x'
	v = ones(length(x))
	
	L = -sparse_log(P)
	L2 = P.*L.^2
	L = P.*L
	vx = v*x
	S = (I-vx)+(P-vx)*inv(I-P+vx)
	covariance = (x*L*S*L*v)[1]
	variance = (x*L2*v)[1] -(x*L*v)[1]^2
	lyapunov =  variance + 2*covariance
	if imag(covariance) < 1.0e-3
	   return lyapunov, variance, covariance, :Success
	else
	   @show covariance
	   return real(lyapunov), real(variance), real(covariance), :ComplexCovariancveWarning
	end
 end
 
