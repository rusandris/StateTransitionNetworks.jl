"""
	timeseries_to_grid(timeseries, grid) -> cell_coordinates,vertex_names
Discretizes a 2D timeseries/trajectory on a grid. Returns
a discrete timeseries containing the the cell coordinates and the list
of vertices with cell coordinates. 
"""
function timeseries_to_grid(timeseries, grid)    
    M = zeros(grid,grid)
    T = length(timeseries[:,1])
    x_min = minimum(timeseries[:, 1])
    y_min = minimum(timeseries[:, 2])
    x_max = maximum(timeseries[:, 1])
    y_max = maximum(timeseries[:, 2])
    dx = 0.5*(x_max-x_min)/grid
    dy = 0.5*(y_max-y_min)/grid
    x_grid = range(x_min-dx, x_max+dx, grid);
    x_min = x_grid[1]
    y_grid = range(y_min-dy, y_max+dy, grid);
    y_min = y_grid[1]
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


"""
	create_stn(discrete_timeseries, vertex_names) -> stn
Creates a state transition network (STN) using the discrete timeseries and vertex list.
The network is a directed metagraph object
## Vertex properties
	stn[i] -> position::Tuple{Int64, Int64}, 
## Edge properties 
	stn[i,j] -> (prob,weigth)::Tuple{Float64,Float64}
"""
function create_stn(discrete_timeseries,vertex_names)
	nr_vertices = length(vertex_names[:,1])
	states = zeros(Int32,length(discrete_timeseries[:,1])) #discrete timeseries but only with vertex indices
	
	for v in 1:nr_vertices
		timepoints_at_v = findall(x -> x == vertex_names[v,2:end],[eachrow(discrete_timeseries)...]) #find all times when state was v
		states[timepoints_at_v] .= v  
	end		
	
	#weight and transition probability matrices
	Q = zeros(nr_vertices, nr_vertices)
    P = zeros(nr_vertices, nr_vertices)
    
    #count transitions
    next_states = circshift(states,-1)
    for i in eachindex(states[1:end-1])
        Q[states[i],next_states[i]] += 1
    end

	#normalize Q and fill P by normalizing rows
    Q = Q./sum(Q)
    for i in 1:size(Q)[1]
        P[i,:] = Q[i,:]./sum(Q[i,:])
    end

	#create directed metagraph with static label and metadata types and default weight 0
	stn = MetaGraph(DiGraph(),Label = Int64, 
		VertexData = Tuple{Int64, Int64}, 
		EdgeData = Tuple{Float64,Float64}, 
		default_weight = 0.0)

	
	#add edges and properties
	#Properties: vertices -> position::Tuple{Int64, Int64}, edges -> (weigths,prob)::Tuple{Float64,Float64}
	
	for v in 1:nr_vertices
		position = tuple(vertex_names[v,2:end]...) 
		stn[v] = position
	end
	
	for i in 1:nr_vertices
        for j in 1:nr_vertices
            if Q[i, j] != 0
				stn[i,j] = (P[i,j],Q[i,j])
            end
        end
    end
    
	return stn
end

