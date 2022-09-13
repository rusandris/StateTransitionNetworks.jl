function timeseries_to_grid(timeseries, grid)
    """
    Discretizes a 2D timeseries/trajectory on a grid. Returns
    a discrete timeseries containing the the cell coordinates and the list
    of vertices with cell coordinates. 
    """
    
    M = zeros(grid,grid)
    T = length(timeseries[:,1])
    x_min = minimum(timeseries[:, 1])
    y_min = minimum(timeseries[:, 2])
    x_max = maximum(timeseries[:, 1])
    y_max = maximum(timeseries[:, 2])
    x_grid = range(x_min, x_max, grid);
    y_grid = range(y_min, y_max, grid);
    x_n = Vector{Int64}(undef, T)
    y_n = Vector{Int64}(undef, T)
    num_vertex = 0
    vertex_names = [];
    x_n, y_n = [], [];

    for row in eachrow(timeseries)
        y = floor(Int,(row[2]-y_min)/Float64(y_grid.step))+1
        x = floor(Int,(row[1]-x_min)/Float64(x_grid.step))+1
        if M[y,x] == 0
            num_vertex += 1 
            push!(vertex_names, [num_vertex, x, y])
            M[y,x] = 1
        end
        push!(x_n, x)
        push!(y_n, y)
    end
    vertex_names = reduce(hcat, vertex_names)'
    d_timeseries = [x_n y_n]
    return d_timeseries, vertex_names
end


function create_STN(discrete_timeseries, vertex_names)
    """
    Creates a state transition network (STN) using the discrete timeseries and vertex list.
    The network is a SimpleWeightedDiGraph object. Two graphs are returned,
    one with occurence probability (Q_ij) and one with the transition probability
    as weights (P_ij) 
    """
    steps = []

    for row in eachrow(discrete_timeseries)
        index = Int[ row == vertex_names[i,2:end] for i in eachindex(vertex_names[:,1])]
        step = vertex_names[findfirst(x -> x == 1, index),1]  
        push!(steps,step)
    end

    num_vertex = length(vertex_names[:,1])

    Q = zeros(num_vertex, num_vertex);
    P = zeros(num_vertex, num_vertex);
    next_step = circshift(steps,-1)

    for i in eachindex(steps[1:end-1])
        Q[steps[i],next_step[i]] += 1
    end

    Q = Q./sum(Q)
    for i in 1:size(Q)[1]
        P[i,:] = Q[i,:]./sum(Q[i,:])
    end

    stn_weight = SimpleWeightedDiGraph(num_vertex)
    stn_prob = SimpleWeightedDiGraph(num_vertex)

    for i in 1:num_vertex
        for j in 1:num_vertex
            if Q[i, j] != 0
                add_edge!(stn_weight, i, j, Q[i,j])
                add_edge!(stn_prob, i, j, P[i, j])
            end
        end
    end

    return stn_weight, stn_prob
end

