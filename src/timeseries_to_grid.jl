export timeseries_to_grid,timeseries_to_grid!,assign_grid_cell
TimeSeries{T,D} = Union{StateSpaceSet{D,T},Array{T,D}} where {T,D}

#--------------------------------------higher level timeseries_to_grid methods (work for 1D and 2D)--------------------------------------

function timeseries_to_grid!(symbolic_timeseries::Vector{I},timeseries::TimeSeries, grid_size::Int64; grid_edges::Vector{Float64} = Float64[],nfloat::Int=100,outside_grid::Symbol=:error) where I<:Integer    

	#decide if grid_edges are needed 
	if isempty(grid_edges)
		#max and min are calculated, grid_edges are determined
		grid_steps,grid_edges = determine_grid_step(timeseries,grid_size,nfloat)
	else
		grid_steps = determine_grid_step(grid_edges,grid_size)
	end

	#select method for handling points falling outside the grid
	if outside_grid == :include
		method = include_outside_grid
	elseif outside_grid == :skip
		throw(ArgumentError("Skipping is not allowed when using in-place methods. Use out-of-place method instead."))
	else
		method = error_outside_grid
	end

	#here gets dispatched to 1D or 2D
	timeseries_to_grid!(symbolic_timeseries,timeseries,grid_size,grid_edges,grid_steps;method=method)

end

function timeseries_to_grid(timeseries::TimeSeries, grid_size::Int64; grid_edges::Vector{Float64} = Float64[],nfloat::Int=100,outside_grid::Symbol=:error)    

	#decide if grid_edges are needed 
	if isempty(grid_edges)
		#max and min are calculated, grid_edges are determined
		grid_steps,grid_edges = determine_grid_step(timeseries,grid_size,nfloat)
	else
		grid_steps = determine_grid_step(grid_edges,grid_size)
	end

	#select method for handling points falling outside the grid
	if outside_grid == :include
		method = include_outside_grid
	elseif outside_grid == :skip
		method = skip_outside_grid
	else
		method = error_outside_grid
	end
	#here gets dispatched to 1D or 2D
	return timeseries_to_grid(timeseries,grid_size,grid_edges,grid_steps;method=method)
end


#------------------------------------------1D methods----------------------------------------

#partitioning into regular equally sized bins

function timeseries_to_grid(timeseries::TimeSeries{T,1}, grid_size::Int,grid_edges::Vector{Float64},grid_steps::Vector{Float64};method=error_outside_grid) where T <:AbstractFloat
   
	if length(grid_edges) != 2
		throw(ArgumentError("Length of `grid_edges` should be 2 in case of 1D time series, got $(length(grid_edges)) !"))
	end

	L, = size(timeseries)

	if method == skip_outside_grid
		symbolic_timeseries = Int[]
	else
		symbolic_timeseries = Vector{Int64}(undef, L)
	end
	
	for p in 1:L
		cell = method(timeseries[p][1],grid_edges,grid_size,grid_steps[1])
		if method == skip_outside_grid
			if cell == 0 
				continue
			end
			push!(symbolic_timeseries,cell)
		else
			symbolic_timeseries[p] = cell
		end
	end
	return symbolic_timeseries
end

function timeseries_to_grid!(symbolic_timeseries::Vector{I},timeseries::TimeSeries{T,1}, grid_size::Integer,grid_edges::Vector{Float64},grid_steps::Vector{Float64};method=error_outside_grid)  where {I<:Integer,T<:AbstractFloat}

	if length(grid_edges) != 2
		throw(ArgumentError("Length of `grid_edges` should be 2 in case of 1D time series, got $(length(grid_edges)) !"))
	end

	min,max = grid_edges
	L, = size(timeseries) 

	for p in 1:L
		cell = method(timeseries[p][1],[min,max],grid_size,grid_steps[1])
		symbolic_timeseries[p] = cell
	end
end


function assign_grid_cell(x::T,x_min::Float64,grid_step::Float64) where T <:AbstractFloat
	floor(Int,(x-x_min)/grid_step) + 1
end


function minmaxima(ts::Array{T,N}) where {T <: AbstractFloat, N}  
	extremas = extrema(ts;dims=1)
	mins = zeros(length(extremas))
	maxs = zeros(length(extremas))

	for i in 1:length(extremas)
		mins[i] = extremas[i][1]
		maxs[i] = extremas[i][2]
	end
	return mins,maxs
end

#-------------------------------more general partitioning (experimental)----------------------------------------

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

#------------------------------determine_grid_step methods--------------------------


function determine_grid_step(grid_edges::Vector{Float64},grid_size::Int)
	dim = Int(length(grid_edges)/2)

	if dim != 1 && dim != 2
		throw(ArgumentError("Length of grid edges list isn't correct (1D and 2D are allowed)!"))
	end

	grid_steps = zeros(dim)

	for i in 1:dim 
		max = grid_edges[i+dim] #pick out max of the grid along a given axis
		min = grid_edges[i] #pick out min of the grid along a given axis
		#grid cell width
		grid_steps[i] = (max - min)/grid_size
	end
	return grid_steps
end


function determine_grid_step(timeseries::TimeSeries,grid_size::Int,nfloat::Int=100)
	#determine mins and maxs
	mins,maxs = minmaxima(timeseries)
	grid_edges = [mins...,maxs...]

	grid_border_iter = 1
	dim = Int(length(grid_edges)/2)
	grid_steps = zeros(dim)

	maxs_plus = zeros(dim)
	max_symbols = zeros(Int,dim)

	#iterate until grid border padding is enough
	while true
		#iterate through axes
		for i in 1:dim 
			max = grid_edges[i+dim] #pick out max of the grid along a given axis
			min = grid_edges[i] #pick out min of the grid along a given axis

			#add padding to right border of the grid
			#so the points on the edge aren't singled out 
			maxs_plus[i] = nextfloat(max,grid_border_iter*nfloat)
			#grid cell width
			grid_steps[i] = (maxs_plus[i] - min)/grid_size
			#test if largest symbol goes out of bounds
			max_symbols[i] = assign_grid_cell(max,min,grid_steps[i])
		end

		if all(s -> s <= grid_size,max_symbols)
			break
		end
		grid_border_iter+=1

	end
	return grid_steps,grid_edges
end

#-------------------------------------------2D timeseries_to_grid methods----------------------------------------

function timeseries_to_grid!(symbolic_timeseries,timeseries::TimeSeries{T,2},grid_size::Int,grid_edges::Vector{Float64},grid_steps::Vector{Float64};method=error_outside_grid) where T <:AbstractFloat   

	if length(grid_edges) != 4
		throw(ArgumentError("Length of `grid_edges` should be 4 in case of 2D time series, got $(length(grid_edges)) !"))
	end

	x_min,y_min,x_max,y_max = grid_edges
	L, = size(timeseries) 

	for p in 1:L
		#handle points outside on the x-axis
		j_cell = method(timeseries[p,1],[x_min,x_max],grid_size,grid_steps[1])
		#handle points outside on the y-axis
		i_cell = method(timeseries[p,2],[y_min,y_max],grid_size,grid_steps[2])
		symbolic_timeseries[p] = (j_cell-1)*grid_size + i_cell 
	end

end

function timeseries_to_grid(timeseries::TimeSeries{T,2},grid_size::Int,grid_edges::Vector{Float64},grid_steps::Vector{Float64};method=error_outside_grid)  where T <:AbstractFloat 

	if length(grid_edges) != 4
		throw(ArgumentError("Length of `grid_edges` should be 4 in case of 2D time series, got $(length(grid_edges)) !"))
	end

	x_min,y_min,x_max,y_max = grid_edges
	L, = size(timeseries)
	
	if method == skip_outside_grid
		symbolic_timeseries = Int[]
	else
		symbolic_timeseries = Vector{Int64}(undef, L)
	end
	
	for p in 1:L
		#handle points outside on the x-axis
		j_cell = method(timeseries[p,1],[x_min,x_max],grid_size,grid_steps[1])
		#handle points outside on the y-axis
		i_cell = method(timeseries[p,2],[y_min,y_max],grid_size,grid_steps[2])

		if method == skip_outside_grid
			if j_cell == 0 || i_cell == 0
				continue
			end
	
			symbol = (j_cell-1)*grid_size + i_cell 
			push!(symbolic_timeseries,symbol)
		else
			symbolic_timeseries[p] = (j_cell-1)*grid_size + i_cell 
		end

	end
	return symbolic_timeseries
end

#--------------------------------------handling points outside of the grid--------------------------------------------

#-------------------------------error_outside_grid--------------------------

function error_outside_grid(v,edges,grid_size,grid_step)
	min,max =  edges 

	if is_inside(v,min,max)
		return assign_grid_cell(v,min,grid_step)
	else
		throw(DomainError(v, "Value falls outside of the grid partition edges " * "[" *" $min , $max " * "]" ))
	end
end

#-------------------------------include_outside_grid--------------------------

function include_outside_grid(v,edges,grid_size,grid_step)
	min,max =  edges 
	if v < min
		cell_index = 1 	
	elseif v > max 
		cell_index = grid_size
	else
		cell_index = assign_grid_cell(v,min,grid_step)
	end
	return cell_index
end

#-------------------------------skip_outside_grid--------------------------

function skip_outside_grid(v,edges,grid_size,grid_step)
	min,max =  edges 

	if is_inside(v,min,max)
		return assign_grid_cell(v,min,grid_step)
	else
		#use 0 symbol to signal values outside the grid
		return 0
	end
end

function is_inside(x,min,max) 
	return min <= x <= max
end
