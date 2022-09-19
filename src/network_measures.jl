function random_walk_on_weighted_graph(graph, N_steps, transient)
    """
    Conducts a random walk process on a weighted graph for `N_steps` discarding
    the first `transient` steps. Returns the normalized walk length.
    """
    n = 0;
    source = sample(1:nv(graph));
    walk_length = 0.0;

    while n < N_steps
        neigh = outneighbors(graph, source)
        w = Vector{Float64}(undef, length(neigh))
        for i in eachindex(neigh)
            w[i] = get_weight(graph, source, neigh[i])
        end
        source = sample(neigh, Weights(w))
        if n > transient
            l = -log(w[indexin(source, neigh)][1])
            walk_length = walk_length + l
        end
        n += 1
    end
    return walk_length
end

function walk_statistics(ensemble, graph, N_steps, transient=5000)
    """
    Calculates the Sinai-Kolmogorov Entropy and Lyapunov measure of a STN
    by calculating the average walk length and the variance of walk lenghts
    over an ensemble of random walks on weighted graph
    """
    walk_length = Vector{Float64}(undef, ensemble)
    for i in 1:ensemble
        walk_length[i] = random_walk_on_weighted_graph(graph, N_steps, transient)
    end
   	entropy = mean(walk_length)/(N_steps-transient)
    lyapunov_measure = var(walk_length)/(N_steps-transient)
    return entropy, lyapunov_measure
end

"""
Calculates analytically the Sinai-Kolmogorov entropy on a STN. The main input `graph_P` is a SimpleWeightedDiGraph
object containg the occurence probability of each transition (Q_ij). The additional `graph_P` argument is another 
graph containg the normalized transition probabilities (P_ij). If this is not provided the function calculates automatically these transition
probabilities.
"""
function sinai_kolmogorov_entropy(graph_Q; graph_P=nothing)
    Q = Matrix(weights(graph_Q))
    if isnothing(graph_P)
        P = Matrix(weights(graph_Q))
        for i in 1:size(Q)[1]
            P[i,:]./=sum(Q[i,:])
        end
    else
        P = Matrix(weights(graph_P))
    end
    entropy = -sum(Q[Q .!=0] .* log.(P[P .!=0]))
    return entropy
end