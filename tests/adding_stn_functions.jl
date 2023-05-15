import StateTransitionNetworks.lyapunov_measure
import SparseArrays: spzeros 
import StateTransitionNetworks.renormalize
import Graphs: DiGraph


function timeseries_to_common_grid(timeseries, grid, x_min, x_max, y_min, y_max)    
    T = length(timeseries[:,1])
    dx = 0.5*(x_max-x_min)/grid
    dy = 0.5*(y_max-y_min)/grid
    x_grid = range(x_min-dx, x_max+dx, grid);
    x_min = x_grid[1]
    y_grid = range(y_min-dy, y_max+dy, grid);
    y_min = y_grid[1]
    x_n = Vector{Int64}(undef, T)
    y_n = Vector{Int64}(undef, T)
    x_n, y_n = [], [];

    for row in eachrow(timeseries)
        y = floor(Int,(row[2]-y_min)/Float64(y_grid.step))+1
        x = floor(Int,(row[1]-x_min)/Float64(x_grid.step))+1
        push!(x_n, x)
        push!(y_n, y)
    end
    d_timeseries = [x_n y_n]
    return d_timeseries
end



function add_timeseries(dt_list, grid; make_ergodic=false, verbose=false)
    Q_null = zeros(Int32, grid*grid, grid*grid) # Null matrix with all possible transitions
    vertex_names = [];  # Future name of vertices
    vertex_place = [];  # Rows, and column number in Q_null for a given vertex
    M = zeros(grid,grid)
    nr_vertices = 0;
    for dt in dt_list
        states = []
        O = zeros(Int32, grid*grid, grid*grid)
        for row in eachrow(dt)
            v = (row[1]-1)*grid + row[2]
            x, y = row
            if M[y,x] == 0
                nr_vertices += 1 
                push!(vertex_names, [nr_vertices, x, y])
                push!(vertex_place, v)
                M[y,x] = 1
            end    
            push!(states, v) 
        end
        #count transitions
        next_states = circshift(states,-1)
        for i in eachindex(states[1:end-1])
            O[states[i],next_states[i]] += 1
        end
        Q_null = Q_null + O
    end
    vertex_names = reduce(hcat, vertex_names)'
    Q = Q_null[vertex_place, vertex_place]
	#weight and transition probability matrices
    P = spzeros(nr_vertices, nr_vertices)
	#normalize Q and fill P by normalizing rows
    Q = Q./sum(Q)

	P = renormalize(Q)
	#create directed metagraph with static label and metadata types and default weight 0
	stn, ret_code = create_stn(P; make_ergodic=make_ergodic, verbose=verbose)
	
    for v in values(stn.vertex_labels)
		x,y = vertex_names[v,2:end] 
		stn[v] = Dict(:x => x,:y => y)
	end
	
    return stn, ret_code
end

function get_grid_edges(psections)
    x_max, y_max = maximum(reduce(vcat, maximum.(Matrix.(psections); dims=1)); dims=1)
    x_min, y_min = minimum(reduce(vcat, minimum.(Matrix.(psections); dims=1)); dims=1)
    return x_min, x_max, y_min, y_max
end

