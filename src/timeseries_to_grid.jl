export timeseries_to_grid,timeseries_to_grid!
TimeSeries{T,D} = Union{StateSpaceSet{D,T},AbstractArray{T,D}} where {T,D}

"""
	timeseries_to_grid(timeseries, grid_size) -> symbolic_timeseries
Discretizes a 2D timeseries/trajectory on a grid.
"""
function timeseries_to_grid(timeseries::TimeSeries{T,2}, grid_size::Integer; grid_edges = [],return_vertex_positions=false) where T <: Real   
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
	
	if return_vertex_positions
		vertex_positions = Dict{Int64,Tuple{Int64, Int64}}()
	end

	#would be nice to write a function for the symbol conversion formula

    for (i,row) in enumerate(timeseries)
        y = floor(Int,(row[2]-y_min)/Float64(y_grid_step))
        x = floor(Int,(row[1]-x_min)/Float64(x_grid_step))
		symbolic_timeseries[i] = x*grid_size + y + 1
		
		if return_vertex_positions
			push!(vertex_positions, symbolic_timeseries[i] => (x,y) )
		end
		
    end
    
    if return_vertex_positions
	    return symbolic_timeseries,vertex_positions
	end
    
    return symbolic_timeseries
end

function timeseries_to_grid(timeseries::TimeSeries{T,1}, grid_size::Integer; grid_edges = [],return_vertex_positions=false)  where T <: Real    
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
    
    if return_vertex_positions
    	vertex_positions = Dict{Int64,Tuple{Int64, Int64}}([symbol => (symbol,0) for symbol in symbolic_timeseries])
    	return symbolic_timeseries,vertex_positions
    end
    return symbolic_timeseries
end

function timeseries_to_grid!(symbolic_timeseries::TimeSeries{T1,1},timeseries::TimeSeries{T2,1}, grid_size::Integer; grid_edges::AbstractVector = Float64[])  where T1 <: Integer where T2 <: Real
    L = length(timeseries)

    if isempty(grid_edges)		
		x_min,x_max = extrema(timeseries)
    else
    	x_min,x_max = grid_edges 
    end
	
	x_max_plus = nextfloat(x_max,100)
   
	x_grid_step = (x_max_plus - x_min)/grid_size

    for i in 1:L
		symbolic_timeseries[i] = floor(Int,(timeseries[i][1]-x_min)/x_grid_step) + 1
    end
	return nothing
end

function extrema(s::StateSpaceSet)
	min,max = minmaxima(s)
	return min[1],max[1]
end

