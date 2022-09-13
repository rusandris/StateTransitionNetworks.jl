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
    walk_length = walk_length./(N_steps-transient)
    entropy = mean(walk_length)
    lyapunov_measure = var(walk_length)
    return entropy, lyapunov_measure
end
