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
	stn[i] -> (:x => x,:y => y)
## Edge properties 
	stn[i,j] -> (:prob => P[i,j],:weight => Q[i,j])
"""
function create_stn(discrete_timeseries,vertex_names)
	nr_vertices = length(vertex_names[:,1])
	states = zeros(Int32,length(discrete_timeseries[:,1])) #discrete timeseries but only with vertex indices
	
	for v in 1:nr_vertices
		timepoints_at_v = findall(x -> x == vertex_names[v,2:end],[eachrow(discrete_timeseries)...]) #find all times when state was v
		states[timepoints_at_v] .= v  
	end		
	
	#weight and transition probability matrices
	Q = spzeros(nr_vertices, nr_vertices)
    P = spzeros(nr_vertices, nr_vertices)
    
    #count transitions
    next_states = circshift(states,-1)
    for i in eachindex(states[1:end-1])
        Q[states[i],next_states[i]] += 1
    end

	#normalize Q and fill P by normalizing rows
    Q = Q./sum(Q)
    for i in 1:size(Q)[1]
        P[i,:] = Q[i,:]./sum(Q[i,:])
       #all(i -> isfinite(i),P[i,:]) || @warn "Inf/NaN values in the transition matrix P[i,j]!"
    end

	#create directed metagraph with static label and metadata types and default weight 0
	stn = @suppress MetaGraph(DiGraph(),Label = Int64, 
		VertexData = Dict{Symbol, Int64}, 
		EdgeData = Dict{Symbol,Float64}, 
		default_weight = 0.0)

	
	#add edges and properties
	#Properties: vertices ->  Dict{Symbol,Int64}, edges -> Dict{Symbol,Float64}
	
	for v in 1:nr_vertices
		x,y = vertex_names[v,2:end] 
		stn[v] = Dict(:x => x,:y => y)
	end
	
	for i in 1:nr_vertices
        for j in 1:nr_vertices
            if P[i, j] != 0
				stn[i,j] = Dict(:prob => P[i,j],:weight => Q[i,j])
            end
        end
    end
    
    #is_strongly_connected(stn) || @warn "The graph is not strongly connected. Increase the length of your timeseries!"
	
	retcode = check_stn(Q,P)
	
	return stn,retcode
end

"""
	create_stn(P) -> stn
Creates a state transition network (STN) using a transition probability matrix.
The network is a directed metagraph object. No error checking in this case.
## Edge properties 
	stn[i,j] -> (:prob => P[i,j])
"""
function create_stn(P)
	nr_vertices = size(P)[1]

	#create directed metagraph with static label and metadata types and default weight 0
	stn = @suppress MetaGraph(DiGraph(),Label = Int64, 
		EdgeData = Dict{Symbol,Float64}, 
		default_weight = 0.0)

	for v in 1:nr_vertices
		stn[v] = nothing
	end
	
	for i in 1:nr_vertices
        for j in 1:nr_vertices
            if P[i, j] != 0
				stn[i,j] = Dict(:prob => P[i,j])
            end
        end
    end
    
    #is_strongly_connected(stn) || @warn "The graph is not strongly connected. Increase the length of your timeseries!"
	
	#retcode = check_stn(Q,P)
	
	return stn	
end

function prob_matrix(stn)
	nr_vertices = nv(stn)
	P = spzeros(Float64,(nr_vertices,nr_vertices))
	for edge in collect(edges(stn))
		i = edge.src
		j = edge.dst
		ilabel = label_for(stn,i)
		jlabel = label_for(stn,j)
		P[i,j] = stn[ilabel,jlabel][:prob]
	end
	return P
end

function weight_matrix(stn)
	nr_vertices = nv(stn)
	Q = spzeros(Float64,(nr_vertices,nr_vertices))
	for edge in collect(edges(stn))
		i = edge.src
		j = edge.dst
		ilabel = label_for(stn,i)
		jlabel = label_for(stn,j)
		Q[i,j] = stn[ilabel,jlabel][:weight]
	end
	return Q
end

function check_stn(Q,P)
	nr_vertices = size(Q)[1]
	default_retcode =:Success
	for i in 1:nr_vertices
		if	sum(Q[:,i]) == 0 
			@warn "Vertex with no incoming edge detected! Length of transient/timeseries is insufficient!"
			return :NoIncoming
		elseif sum(Q[i,:]) == 0 
			@warn "Vertex with no outgoing edge detected! Length of timeseries is insufficient!"
			return :NoOutgoing
		elseif  P[i,i] == 1
			@warn "Dead-end (self-loop edge) detected!"
			return :DeadEnd
		end
	end
	return default_retcode
end



