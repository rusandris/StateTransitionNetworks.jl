"""
	timeseries_to_grid(timeseries, grid) -> cell_coordinates,vertex_names
Discretizes a 2D timeseries/trajectory on a grid. Returns
a discrete timeseries containing the the cell coordinates and the list
of vertices with cell coordinates. 
"""
function timeseries_to_grid(timeseries, grid; grid_edges = [])    
    M = zeros(grid,grid)
    T = length(timeseries[:,1])
    
    if isempty(grid_edges)		
		x_min = minimum(timeseries[:, 1])
		y_min = minimum(timeseries[:, 2])
		x_max = maximum(timeseries[:, 1])
		y_max = maximum(timeseries[:, 2])
    else
    	x_min,y_min,x_max,y_max = grid_edges 
    end
    #partitioning with little extra space on both ends
    
    dx = 0.5*(x_max-x_min)/grid
    dy = 0.5*(y_max-y_min)/grid
    
    x_grid = range(x_min-dx, x_max+dx, grid);
    x_min = x_grid[1]
    y_grid = range(y_min-dy, y_max+dy, grid);
    y_min = y_grid[1]
    
    #arrays for space-discrete timeseries 
    x_n = Vector{Int64}(undef, T)
    y_n = Vector{Int64}(undef, T)
    num_vertex = 0
    vertex_names = OrderedDict{Tuple{Int64,Int64},Int64}()
    x_n, y_n = [], [];

    for row in eachrow(timeseries)
        y = floor(Int,(row[2]-y_min)/Float64(y_grid.step))+1
        x = floor(Int,(row[1]-x_min)/Float64(x_grid.step))+1
        if M[y,x] == 0
            num_vertex += 1 
            vertex_names[(x,y)] = num_vertex
            M[y,x] = 1
        end
        push!(x_n, x)
        push!(y_n, y)
    end
    d_timeseries = collect(zip(x_n,y_n))
    return d_timeseries, vertex_names
end

"""
	create_stn(ts,grid::Int64,plane,idxs;make_ergodic=false, verbose=false,kwargs...) -> stn
Higher level function that creates a state transition network (STN) directly from the timeseries of a continuous system `ts` using `ChaosTools.poincaresos`. The returned `stn` is a directed metagraph object. \\
`plane` is the same as in `ChaosTools.poincaresos`. With `idxs` you can choose the variables you want to save.
`make_ergodic=true` returns an STN with only one strongly connected components. Defaults to `false`.   
`verbose=true` logs info about the network checking process. Defaults to `false`.
Other `kwargs` get propageted into `ChaosTools.poincaresos`.

For more info about PSOS  go to https://juliadynamics.github.io/DynamicalSystems.jl/v1.3/chaos/orbitdiagram/#ChaosTools.poincaresos
"""
function create_stn(ts,grid::Int64,plane,idxs;make_ergodic=false, verbose=false,kwargs...)
	
	#psos
	psection = poincaresos(Dataset(ts), plane; idxs=idxs,warning=true,kwargs...)
	#method for time-discrete trajectory
	create_stn(psection,grid; make_ergodic=make_ergodic, verbose=verbose)
	
end


"""
	create_stn(time_discrete_ts, vertex_names;make_ergodic=false,verbose=false) -> stn
Higher level function that creates a state transition network (STN) directly from the  time-discrete timeseries (could be a trajectory of a discrete dynamical system). The returned `stn` is a directed metagraph object. \\
`make_ergodic=true` returns an STN with only one strongly connected components. Defaults to `false`.   
`verbose=true` logs info about the network checking process. Defaults to `false`.  

For more info about the network checking go to Graphs.jl: https://juliagraphs.org/Graphs.jl/dev/algorithms/connectivity/#Graphs.strongly_connected_components 
"""
function create_stn(time_discrete_ts,grid::Int64; make_ergodic=false, verbose=false)
	discrete_timeseries,vertex_names = timeseries_to_grid(time_discrete_ts,grid)
	create_stn(discrete_timeseries,vertex_names;make_ergodic=make_ergodic,verbose=verbose)
end

"""
	create_stn(discrete_timeseries, vertex_names;make_ergodic=false,verbose=false) -> stn
Creates a state transition network (STN) using the discrete timeseries and vertex list. The network is a directed metagraph object. \\
`make_ergodic=true` returns an STN with no defective vertices (vertices with no ingoing/outgoing edges and deadends are deleted). Defaults to `false`.  

## Vertex properties
	stn[i] -> (:x => x,:y => y)
## Edge properties 
	stn[i,j] -> (:prob => P[i,j],:weight => Q[i,j])
"""
function create_stn(discrete_timeseries,vertex_names;make_ergodic=false,verbose=false)
	nr_vertices = length(vertex_names)
	
	#weight and transition probability matrices
	Q = spzeros(nr_vertices, nr_vertices)
    P = spzeros(nr_vertices, nr_vertices)
    
    #probability distribution of states
    states_distrib = zeros(nr_vertices)
    
    #count transitions
    for i in eachindex(discrete_timeseries[1:end-1])
    	state = discrete_timeseries[i]
    	next_state = discrete_timeseries[i+1]
        Q[vertex_names[state],vertex_names[next_state]] += 1
        states_distrib[vertex_names[state]] += 1
    end
    
    #update end (final) state 
    end_state = discrete_timeseries[end]
    states_distrib[vertex_names[end_state]] += 1
    
    #normalize state distribution
	states_distrib = states_distrib ./ sum(states_distrib)

	#normalize Q and fill P by normalizing rows
    Q = Q./sum(Q)
	P = renormalize(Q)

	#create directed metagraph with static label and metadata types and default weight 0
	stn = MetaGraph(
        DiGraph(),
        Int64,
        Dict{Symbol, Union{Int64,Float64}},
        Dict{Symbol, Float64},
        nothing,
        edge_data -> 1.0,
        0.0)

	
	#add edges and properties
	#Properties: vertices ->  Dict{Symbol,Int64}, edges -> Dict{Symbol,Float64}
	
	for state in keys(vertex_names)
		stn[vertex_names[state]] = Dict(:x => state[1],:y => state[2],:prob => states_distrib[vertex_names[state]])
	end
	
	for i in 1:length(vertex_names)
        for j in 1:length(vertex_names)
            if P[i, j] != 0
				stn[i,j] = Dict(:prob => P[i,j],:weight => Q[i,j])
            end
        end
    end

	retcode = check_stn!(stn,make_ergodic=make_ergodic,verbose=verbose)
	
	return stn,retcode
end

"""
	create_stn(P) -> stn
Creates a state transition network (STN) using a transition probability matrix.
The network is a directed metagraph object. No error checking in this case.
## Edge properties 
	stn[i,j] -> (:prob => P[i,j])
"""
function create_stn(P::AbstractMatrix;make_ergodic=false,verbose=false)
	renormalize!(P)

	nr_vertices = size(P)[1]

	#create directed metagraph with static label and metadata types and default weight 0
	stn = MetaGraph(
        DiGraph(),
        Int64,
        Dict{Symbol, Int64},
        Dict{Symbol, Float64},
        nothing,
        edge_data -> 1.0,
        0.0)

	for v in 1:nr_vertices
		stn[v] = Dict{Symbol, Int64}()
	end
	
	Q = calculate_weight_matrix(P)
	
	for i in 1:nr_vertices
        for j in 1:nr_vertices
            if P[i, j] != 0
				stn[i,j] = Dict(:prob => P[i,j], :weight => Q[i,j])
            end
        end
    end
    
	retcode = check_stn!(stn;make_ergodic=make_ergodic,verbose=verbose)
	
	return stn,retcode	
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

function renormalize!(stn)
	nr_vertices = nv(stn)
	P = prob_matrix(stn)
	renormalize!(P)
	Q = calculate_weight_matrix(P)
    for i in 1:nr_vertices
        for j in 1:nr_vertices
            if P[i, j] != 0
            	ilabel = label_for(stn,i)
				jlabel = label_for(stn,j)
				stn[ilabel,jlabel] = Dict(:prob => P[i,j],:weight => Q[i,j])
            end
        end
    end

end

function renormalize!(P::AbstractMatrix)
	for i in 1:size(P)[1]
		sumPi = sum(P[i,:])
		if sumPi != 0
    		P[i,:] = P[i,:]./sumPi
    	end
        all(i -> isfinite(i),P[i,:]) || @warn "Stochastic matrix cannot be normalized. Inf/NaN values in the transition matrix P[i,j]!"
    end
end

function renormalize(Q::AbstractMatrix)
	P = spzeros(size(Q))
	for i in 1:size(Q)[1]
		sumQi = sum(Q[i,:])
		if sumQi != 0
    		P[i,:] = Q[i,:]./sumQi
    	else
    		P[i,:] = Q[i,:]
    		#@warn "Stochastic matrix cannot be normalized."
    	end
        all(i -> isfinite(i),P[i,:]) || @warn "Stochastic matrix cannot be normalized. Inf/NaN values in the transition matrix P[i,j]!"
    end
    return P
end

function calculate_weight_matrix(P::AbstractMatrix)
	Q = spzeros(size(P))
	λ, X = eigen(transpose(Matrix(P)))
	if real(λ[end]) ≈ 1  
		x = real(transpose(X[:,end]./sum(X[:,end])))
	else
		@warn "No eigenvalues close to 1. Couldn't calculate the weight matrix."
		return Q  
	end
	
	for i in 1:size(P)[1]
		Q[i,:] = x[i] .* P[i,:]
	end
	Q
end

function isnormalized(P::AbstractMatrix)
	for r in eachrow(P)
		!(sum(r) ≈ 1) && return false
	end
	return true
end


function check_stn!(stn;make_ergodic=false,verbose=false)
	comps = strongly_connected_components(stn) #components list 
	comps_labels = [label_for.(Ref(stn),comp) for comp in comps] #components list with labels from metagraph
	
	nr_comps = length(comps)
	nr_vertices0 = nv(stn)
	if nr_comps == 1
		return :Success
	else
		verbose && @warn "STN is not strongly connected! $nr_comps"*" strongly connected components were found."
		if make_ergodic
			lengths_of_comps = length.(comps)
			largest_component_nrvertices,largest_component_index = findmax(lengths_of_comps)
			
			nr_deleted_vertices = 0 
			for (i,comp) in enumerate(comps_labels)
				if i == largest_component_index
					continue
				else
					for v in comp
						delete!(stn,v)
						nr_deleted_vertices += 1
					end
				end
			end
			remaining_percentage = round((nr_vertices0 - nr_deleted_vertices)/nr_vertices0 * 100; digits=2)
			verbose &&  @info "Using component with $largest_component_nrvertices vertices. Deleted $nr_deleted_vertices"*" out of $nr_vertices0"*" vertices"*" (remaining $remaining_percentage"*"%)"
			if nr_vertices0 - nr_deleted_vertices <= nr_deleted_vertices 
				@warn "Too many deleted vertices! Unusable STN."
				return :Unusable
			end
			
			#rebuild, renormalize here!!!!
			renormalize!(stn)
			return :Success
		else
			return :NotConnected
		end
	end	
end

