

"""
	timeseries_to_grid(timeseries, grid_size) -> cell_coordinates,vertex_names
Discretizes a 2D timeseries/trajectory on a grid. Returns
a discrete timeseries containing the the cell coordinates and the list
of vertices with cell coordinates. 
"""
function timeseries_to_grid(timeseries::TimeSeries, grid_size::Integer; grid_edges = [])    
    M = spzeros(grid_size,grid_size)
    T = length(timeseries[:,1])
    
    if isempty(grid_edges)		
		x_min = minimum(timeseries[:, 1])
		y_min = minimum(timeseries[:, 2])
		x_max = maximum(timeseries[:, 1])
		y_max = maximum(timeseries[:, 2])
    else
    	x_min,y_min,x_max,y_max = grid_edges 
    end
    
	#partitioning between [`x_min`, `x_max`] into `grid` number of cells
    x_grid = range(x_min, nextfloat(x_max, 100), grid_size+1);
    x_min = x_grid[1]
    y_grid = range(y_min, nextfloat(y_max, 100), grid_size+1);
    y_min = y_grid[1]
   
    #arrays for space-discrete timeseries 
    x_n = Vector{Int64}(undef, T)
    y_n = Vector{Int64}(undef, T)
    num_vertex = 0
    vertex_positions = OrderedDict{Tuple{Int64,Int64},Int64}()
    x_n, y_n = [], [];

    for row in timeseries
        y = floor(Int,(row[2]-y_min)/Float64(y_grid.step))+1
        x = floor(Int,(row[1]-x_min)/Float64(x_grid.step))+1
        if M[y,x] == 0
            num_vertex += 1 
            vertex_positions[(x,y)] = num_vertex
            M[y,x] = 1
        end
        push!(x_n, x)
        push!(y_n, y)
    end
    symbolic_timeseries = collect(zip(x_n,y_n))
    return symbolic_timeseries, vertex_positions
end

"""
	create_stn(ts,grid_size::Int64,plane,idxs;make_ergodic=false, verbose=false,kwargs...) -> stn,retcode
Higher level function that creates a state transition network (STN) directly from the timeseries of a continuous system `ts` using `DynamicalSystemsBase.poincaresos`. The returned `stn` is a directed metagraph object. \\

Keyword arguments:
* `plane` is the same as in `DynamicalSystemsBase.poincaresos`. With `idxs` you can choose the variables you want to save.
* `make_ergodic=true` returns an STN with only one strongly connected component. Defaults to `false`.   
* `verbose=true` logs info about the networks' strongly connected components. Defaults to `false`.
Other `kwargs` get propageted into `DynamicalSystemsBase.poincaresos`.

For more info about PSOS  go to https://juliadynamics.github.io/DynamicalSystems.jl/v1.3/chaos/orbitdiagram/#DynamicalSystemsBase.poincaresos
"""
function create_stn(ts::TimeSeries,grid_size::Integer,plane,idxs;make_ergodic=false, verbose=false,grid_edges=[],kwargs...)
	
	#psos
	psection = DynamicalSystemsBase.poincaresos(DynamicalSystemsBase.StateSpaceSet(ts), plane; save_idxs=idxs,warning=true,kwargs...)

	#method for time-discrete trajectory
	create_stn(psection,grid_size; make_ergodic=make_ergodic, verbose=verbose,grid_edges=grid_edges)
	
end


"""
	create_stn(time_discrete_ts, vertex_names;make_ergodic=false,verbose=false) -> stn,retcode
Higher level function that creates a state transition network (STN) directly from the  time-discrete timeseries (could be a trajectory of a discrete dynamical system). The returned `stn` is a directed metagraph object. \\

Keyword arguments:

* `make_ergodic=true` returns an STN with only one strongly connected component. Defaults to `false`.   
* `verbose=true` logs info about the networks' strongly connected components. Defaults to `false`..  

For more info about the network checking go to Graphs.jl: https://juliagraphs.org/Graphs.jl/dev/algorithms/connectivity/#Graphs.strongly_connected_components 
"""
function create_stn(time_discrete_ts::TimeSeries,grid_size::Integer; make_ergodic=false, verbose=false,grid_edges=[])
	symbolic_timeseries,vertex_positions = timeseries_to_grid(time_discrete_ts,grid_size,grid_edges=grid_edges)
	create_stn(symbolic_timeseries,vertex_positions;make_ergodic=make_ergodic,verbose=verbose)
end

"""
	create_stn(symbolic_timeseries, vertex_positions;make_ergodic=false,verbose=false) -> stn,retcode
Creates a state transition network (STN) using the `symbolic timeseries` and the `vertex positions`. The network is a directed metagraph object. \\
Keyword arguments:
* `make_ergodic` : if `true`, only the largest strongly connected component is returned
* `verbose`: Logs information about the networks' strongly connected components

## Vertex properties
	stn[i] -> (:x => x,:y => y)
## Edge properties 
	stn[i,j] -> (:prob => P[i,j],:weight => Q[i,j])
"""
function create_stn(symbolic_timeseries::TimeSeries,vertex_positions::AbstractDict;make_ergodic=false,verbose=false)
	
	P,Q,state_probabilities = calculate_transition_matrix(symbolic_timeseries;symbol_dictionary=vertex_positions,returnQ=true,return_state_distribution=true)
	
	stn, retcode = create_stn(P; make_ergodic=make_ergodic,
		vertex_positions=vertex_positions,
		state_probabilities=state_probabilities,
		Q=Q,
		verbose=verbose)
	
	return stn,retcode
end


"""
	create_stn(P;kwargs...) -> stn,retcode
Creates a state transition network (STN) using a transition probability matrix (stochastic matrix `P`).
The network is a directed metagraph object.
Optional keyword arguments:
* `make_ergodic` : if `true`, only the largest strongly connected component is returned 
* `vertex_positions::AbstractDict{Tuple{Int64, Int64}, Int64}` : dictionary with positions of vertices in the grid. Keys are the grid positions (x::Int64,y::Int64) and values are the vertex label (v::Int64). Dummy coordinates by default `Dict([(v,0) => v for v in 1:size(P)[1]])`.
* `state_probabilites::Vector{Float64}`: vector of the probability of all states
* `Q::AbstractMatrix`: matrix containig the weights (non-conditional probability) of transitions
* `verbose`: Logs information about the networks' strongly connected components

## Vertex properties 
	stn[i,j] -> (:prob,:x,:y)
	
## Edge properties 
	stn[i,j] -> (:prob,:weight)
"""
function create_stn(P::AbstractMatrix;make_ergodic=false,
	vertex_positions::AbstractDict{Tuple{Int64, Int64}, Int64}=Dict([(v,0) => v for v in 1:size(P)[1]]),
	state_probabilities::Vector{Float64}=fill(NaN,size(P)[1]),	
	Q::AbstractMatrix=fill(NaN,(size(P)[1],size(P)[1])),
	verbose=false)

	(!(all(x -> isnan(x), state_probabilities)) && !(sum(state_probabilities) ≈ 1.0)) && throw(ArgumentError("Probabilities of states must sum up to 1!"))
	(!(all(x -> isnan(x), state_probabilities)) && !(sum(Q) ≈ 1.0)) && throw(ArgumentError("Non-conditional probabilities of transitions (weights) must sum up to 1!"))
	
	if !(isnormalized(P))
		renormalize!(P;verbose=verbose)
	end
	
	nr_vertices = size(P)[1]
	
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
	
	for state in keys(vertex_positions)
		stn[vertex_positions[state]] = Dict(:x => state[1],:y => state[2],:prob => state_probabilities[vertex_positions[state]])
	end
	
	for i in 1:length(vertex_positions)
        for j in 1:length(vertex_positions)
            if P[i, j] > 0
				stn[i,j] = Dict(:prob => P[i,j],:weight => Q[i,j])
            end
        end
    end

	retcode = check_stn!(stn,make_ergodic=make_ergodic,verbose=verbose)
	return stn,retcode	
end


"""
	get_state_distribution(stn) -> prob_states,pos_states
Returns the probability distribution of states and the positions of the vertices that correspond to these states.  
"""
function get_state_distribution(stn)
	prob_states = zeros(nv(stn))	#probability of vertices(states)
	pos_states = zeros(Int32,nv(stn),2)	#position of vertices
	
	for v in collect(vertices(stn))
		# v is the code of the vertex
		v_label = label_for(stn,v)
		prob_states[v] = stn[v_label][:prob] 
		pos_states[v,:] .= [stn[v_label][:x],stn[v_label][:y]] 
	end
	
	return prob_states,pos_states
end



function get_transition_matrix(stn)
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

function get_weight_matrix(stn)
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

function renormalize!(stn;verbose=false)
	nr_vertices = nv(stn)
	P = prob_matrix(stn)
	renormalize!(P;verbose=verbose)
	Q = calculate_weight_matrix(P)
    for i in 1:nr_vertices
        for j in 1:nr_vertices
            if P[i, j] > 0
            	ilabel = label_for(stn,i)
				jlabel = label_for(stn,j)
				stn[ilabel,jlabel] = Dict(:prob => P[i,j],:weight => Q[i,j])
            end
        end
    end

end

function renormalize!(P::AbstractMatrix;verbose=false)
	for i in 1:size(P)[1]
		sumPi = sum(P[i,:])
		P[i,:] = P[i,:]./sumPi
    end
    !all(p -> isfinite(p),P) && verbose && @warn "The matrix is not stochastic!";
end


function calculate_weight_matrix(P::AbstractMatrix)
	Q = spzeros(size(P))
	x = stationary_distribution(P)
	
	for i in 1:size(P)[1]
		Q[i,:] = x[i] .* P[i,:]
	end
	return Q
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
		Q = get_weight_matrix(stn)
		if !(sum(Q) ≈ 1)
			renormalize!(stn;verbose=verbose)
		end
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
			else
				# Rebuild, renormalize
				renormalize!(stn;verbose=verbose)
				return :Success
			end
		else
			return :NotConnected
		end
	end	
end

