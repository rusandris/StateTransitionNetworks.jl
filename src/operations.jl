function get_grid_edges(psections)
    x_max, y_max = maximum(reduce(vcat, maximum.(Matrix.(psections); dims=1)); dims=1)
    x_min, y_min = minimum(reduce(vcat, minimum.(Matrix.(psections); dims=1)); dims=1)
    return [x_min, y_min, x_max, y_max]
end

function timeseries_to_common_grid(list_of_timeseries, grid)
    grid_edges = get_grid_edges(list_of_timeseries)
    
    symbolic_ts_list = []
    
    for ts in list_of_timeseries
    	sym_ts,vnames = timeseries_to_grid(ts, grid; grid_edges = grid_edges)
    	push!(symbolic_ts_list,sym_ts)
    end
	return symbolic_ts_list 
end

function add_timeseries(list_of_timeseries,grid,plane;idxs,make_ergodic=false, verbose=false,kwargs...)

	psections = []
	for ts in list_of_timeseries
		psection = poincaresos(DynamicalSystemsBase.StateSpaceSet(ts), plane; save_idxs=idxs,kwargs...)
		push!(psections,psection)
	end
	
	symbolic_ts_list = timeseries_to_common_grid(psections,grid)
	add_timeseries(symbolic_ts_list,grid;verbose=verbose,make_ergodic=make_ergodic) 

end


function add_timeseries(symbolic_ts_list, grid; make_ergodic=false, verbose=false)
    Q_null = zeros(Int32, grid*grid, grid*grid) # Null matrix with all possible transitions
    vertex_names = OrderedDict{Tuple{Int64,Int64},Int64}() 
    vertex_place = [];  # Rows, and column number in Q_null for a given vertex
    M = zeros(grid,grid)
    nr_vertices = 0;
    
    #probability distribution of states
    states_distrib = zeros(nr_vertices)
    
    for ts in symbolic_ts_list
        O = zeros(Int32, grid*grid, grid*grid)
        
        
        for (i,state) in enumerate(ts)
            x, y = state
            v = (x-1)*grid + y
            
            if M[y,x] == 0
                nr_vertices += 1 
                vertex_names[(x,y)] = nr_vertices
                push!(states_distrib,0.0)
                push!(vertex_place, v)
                M[y,x] = 1
            end
            
            if i < length(ts) 
            	x_next,y_next = ts[i+1]
            	v_next = (x_next-1)*grid + y_next
            	O[v,v_next] +=1
            	states_distrib[vertex_names[state]] += 1
            end
            
        end
        
		#update end (final) state 
		end_state = ts[end]
		states_distrib[vertex_names[end_state]] += 1
		
        Q_null = Q_null + O
    end
    
    #normalize state distribution
	states_distrib = states_distrib ./ sum(states_distrib)
    
    Q = Q_null[vertex_place, vertex_place]
	#weight and transition probability matrices
    P = spzeros(nr_vertices, nr_vertices)
	#normalize Q and fill P by normalizing rows
    Q = Q./sum(Q)

	P = calculate_transition_matrix(Q)
	#create directed metagraph with static label and metadata types and default weight 0
	stn, ret_code = create_stn(P; make_ergodic=make_ergodic, verbose=verbose)
	
  	for state in keys(vertex_names)
		stn[vertex_names[state]] = Dict(:x => state[1],:y => state[2],:prob => states_distrib[vertex_names[state]])
	end
	
    return stn, ret_code
end

function are_equal(stn1, stn2)
    Q1 = weight_matrix(stn1)
    Q2 = weight_matrix(stn2)
    if size(Q1) != size(Q2)
        return false
    end
    vertices_2 = Vector(vertices(stn2));
    permutation = [];
    for i in vertices(stn1)
        for j in vertices_2
            if stn1[i] == stn2[j]
                push!(permutation, j)
                deleteat!(vertices_2, findall(x->x==j,vertices_2))
                break
            end
        end
    end
    return Q1 ≈ Q2[permutation, permutation]
end
