export timeseries_to_grid,timeseries_to_grid_1D
TimeSeries = Union{AbstractStateSpaceSet,AbstractArray}

"""
	timeseries_to_grid(timeseries, grid_size) -> symbolic_timeseries
Discretizes a 2D timeseries/trajectory on a grid.
"""
function timeseries_to_grid(timeseries::TimeSeries, grid_size::Integer; grid_edges = [],return_vertex_positions=false)    
    L = length(timeseries[:,1])
    
    
    if isempty(grid_edges)		
		x_min = minimum(timeseries[:, 1])
		y_min = minimum(timeseries[:, 2])
		x_max = maximum(timeseries[:, 1])
		y_max = maximum(timeseries[:, 2])
    else
    	x_min,y_min,x_max,y_max = grid_edges 
    end
	
	x_max_plus = nextfloat(x_max,100)
   
	x_grid_step = (x_max_plus - x_min)/grid_size
	y_grid_step = (x_max_plus - y_min)/grid_size
   
	symbolic_timeseries = Vector{Int64}(undef, L)

    for (i,row) in enumerate(timeseries)
        y = floor(Int,(row[2]-y_min)/Float64(y_grid_step))
        x = floor(Int,(row[1]-x_min)/Float64(x_grid_step))
		symbolic_timeseries[i] = x*grid_size + y + 1
    end
    return symbolic_timeseries
end

function timeseries_to_grid_1D(timeseries::TimeSeries, grid_size::Integer; grid_edges = [],return_vertex_positions=false)    
    L = length(timeseries)

    if isempty(grid_edges)		
		x_min = minimum(timeseries[:, 1])
		x_max = maximum(timeseries[:, 1])
    else
    	x_min,x_max = grid_edges 
    end
	
	x_max_plus = nextfloat(x_max,100)
   
	x_grid_step = (x_max_plus - x_min)/grid_size

	symbolic_timeseries = Vector{Int64}(undef, L)

    for (i,row) in enumerate(timeseries)
        x = floor(Int,(row[1]-x_min)/Float64(x_grid_step))
		symbolic_timeseries[i] = x + 1
    end
    return symbolic_timeseries
end

