export timeseries_to_grid,timeseries_to_grid!,assign_grid_cell
TimeSeries{T,D} = Union{StateSpaceSet{D,T},AbstractArray{T,D}} where {T,D}

#------------------------------------------1D methods----------------------------------------

#partitioning into regular equally sized bins

function timeseries_to_grid(timeseries::TimeSeries{T,1}, grid_size::Int64; grid_edges::Vector{Float64} = Float64[]) where T <:AbstractFloat
   
	L = length(timeseries)

	symbolic_timeseries = Vector{Int64}(undef, L)
	
    timeseries_to_grid!(symbolic_timeseries,timeseries, grid_size; grid_edges=grid_edges)

    return symbolic_timeseries
end


function timeseries_to_grid!(symbolic_timeseries::Vector{I},timeseries::TimeSeries{T,1}, grid_size::Integer; grid_edges::Vector{Float64} = Float64[])  where {I<:Integer,T<:AbstractFloat}
    L = length(timeseries)

    if isempty(grid_edges)		
		x_min,x_max = extrema(timeseries)
    else
    	x_min,x_max = grid_edges 
    end
	
	x_max_plus = nextfloat(x_max,100)
   
	x_grid_step = (x_max_plus - x_min)/grid_size

    for i in 1:L
		symbolic_timeseries[i] = assign_grid_cell(timeseries[i][1],x_min,x_grid_step)
    end
	return nothing
end

function assign_grid_cell(x::T,x_min::Float64,grid_step::Float64) where T <:AbstractFloat
	floor(Int,(x-x_min)/grid_step) + 1
end

function extrema(s::StateSpaceSet)
	min,max = minmaxima(s)
	return min[1],max[1]
end

# more general partitioning 

function timeseries_to_grid(timeseries::TimeSeries{T,1}, partition::Function)  where T <:AbstractFloat 
    L = length(timeseries)

	symbolic_timeseries = Vector{Int64}(undef, L)
	
    timeseries_to_grid!(symbolic_timeseries,timeseries, partition)
    
    return symbolic_timeseries
end

function timeseries_to_grid!(symbolic_timeseries::Vector{I},timeseries::TimeSeries{T,1}, partition::Function) where {I<:Integer,T <:AbstractFloat}
	for i in 1:length(timeseries)
		symbolic_timeseries[i] = partition(timeseries[i][1])
	end
end

#----------------------------------------2D methods-----------------------------------------

"""
	timeseries_to_grid(timeseries, grid_size) -> symbolic_timeseries
Discretizes a 2D timeseries/trajectory on a grid.
"""
function timeseries_to_grid(timeseries::TimeSeries{T,2}, grid_size::Int64; grid_edges::Vector{Float64} = Float64[]) where T <:AbstractFloat   
	L = size(timeseries)[1]

	symbolic_timeseries = Vector{Int64}(undef, L)
	
	timeseries_to_grid!(symbolic_timeseries,timeseries, grid_size; grid_edges = grid_edges)
    
    return symbolic_timeseries
end

function timeseries_to_grid!(symbolic_timeseries::Vector{I},timeseries::TimeSeries{T,2}, grid_size::Int64; grid_edges::Vector{Float64} = Float64[]) where {I<:Integer,T <:AbstractFloat}    
    L = size(timeseries)[1]

	if isempty(grid_edges)		
		x_min,x_max = extrema(timeseries[:, 1])
		y_min,y_max = extrema(timeseries[:, 2])
	else
		x_min,y_min,x_max,y_max = grid_edges 
	end

	x_max_plus = nextfloat(x_max,100)
	y_max_plus = nextfloat(y_max,100)

	x_grid_step = (x_max_plus - x_min)/grid_size
	y_grid_step = (y_max_plus - y_min)/grid_size

	for i in 1:L
		u = timeseries[i,:]
		y_cell = assign_grid_cell(u[1],x_min,x_grid_step)
		x_cell = assign_grid_cell(u[2],y_min,y_grid_step)
		symbolic_timeseries[i] = (x_cell-1)*grid_size + y_cell 
		
	end
end