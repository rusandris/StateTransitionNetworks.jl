"""
	stn_analysis(timeseries::AbstractDataset;grid,plane,ensemble,N_steps) ->  S, Λ
	
Calculates the Sinai-Kolmogorov Entropy and Lyapunov measure of a STN created from `timeseries`. If the retcode is other than `:Success`, 'NaN's are returned.
## Keyword arguments
* `grid` : size of grid used for discretization 
* `plane` : PSOS plane propagated to `poincaresos` from ChaosTools
* `idxs` : choose which variables to save
* `ensemble` : number of individual random walks on the STN
* `N_steps` : number of steps in each random walk
* `return_stn` : returns only the `stn` and the `retcode`
"""
function stn_analysis(timeseries::Dataset;grid,plane,idxs,ensemble=100,N_steps=1000,return_stn=false)
	dim = size(timeseries)[2]
	psection = poincaresos(timeseries, plane;idxs=idxs); 
	
	discrete_timeseries, vertex_names = timeseries_to_grid(psection,grid)
	stn,retcode = create_stn(discrete_timeseries,vertex_names)
	
	if return_stn
		return stn,retcode
	end
	
	if retcode ==:Success
		return network_measures(stn,ensemble,N_steps)
	else
		return NaN,NaN
	end
end

stn_analysis(timeseries::Matrix;grid,plane,idxs,ensemble=100,N_steps=1000,return_stn=false) = stn_analysis(Dataset(timeseries);grid=grid,plane=plane,idxs=idxs,ensemble=ensemble,N_steps=N_steps,return_stn=return_stn)

"""
	read_bin(filename::String,T::DataType,dim) -> Matrix
Reads data from `filename` with elements of type `T` and shape `(:,dim)`. 

"""
function read_bin(filename::String,T::DataType,dim)
	data = reinterpret(T,read(filename))
	return Matrix(reshape(data,(dim,:))')
end

function ndensity(stn)
	nr_vertices = nv(stn)
	ne(stn)/(nr_vertices*(nr_vertices-1))
end



