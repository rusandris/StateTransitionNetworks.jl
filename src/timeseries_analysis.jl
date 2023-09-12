stn_analysis(timeseries::Matrix;grid,plane,idxs,ensemble=100,N_steps=1000,make_ergodic=false, verbose=false,return_stn=false,use_analytic=false,use_stored_distribution=false) = stn_analysis(DynamicalSystemsBase.StateSpaceSet(timeseries);grid=grid,plane=plane,idxs=idxs,ensemble=ensemble,N_steps=N_steps,make_ergodic=make_ergodic, verbose=verbose,return_stn=return_stn,use_analytic=use_analytic,use_stored_distribution=use_stored_distribution)



"""
	stn_analysis(timeseries::DynamicalSystemsBase.StateSpaceSet;grid,plane,idxs,ensemble=100,N_steps=1000,make_ergodic=false, verbose=false,return_stn=false,use_analytic=false) ->  S, Λ
	
Calculates the Sinai-Kolmogorov Entropy and Lyapunov measure of a STN created from `timeseries`. If `retcode` is other than `:Success`, `NaN`s are returned.
## Keyword arguments
* `grid` : size of grid used for discretization 
* `plane` : PSOS plane propagated to `poincaresos` from ChaosTools
* `idxs` : choose which variables to save
* `ensemble` : number of individual random walks on the `stn`
* `N_steps` : number of steps in each random walk
* `return_stn` : returns the `stn` and the `retcode` as well
* `make_ergodic` : returns an `stn` with only one strongly connected component. Defaults to `false`.  
* `verbose` : logs the connectedness checking process
* `use_analytic` : use analytic formula for calculating network measures
* `use_stored_distribution` : use the already stored stationary distribution
"""

function stn_analysis(timeseries::DynamicalSystemsBase.StateSpaceSet;grid,plane,idxs,ensemble=100,N_steps=1000,make_ergodic=false, verbose=false,return_stn=false,use_analytic=false,use_stored_distribution=false)

	stn,retcode = create_stn(timeseries,grid::Int64,plane,idxs;make_ergodic=make_ergodic, verbose=verbose)
	
	entropy,lyapunov = NaN,NaN #initialize variables
	
	if retcode ==:Success
		if use_analytic
			P = prob_matrix(stn)
			entropy,lyapunov = use_stored_distribution ? network_measures(P;x=state_distribution(stn)[1]) : network_measures(P)
		else
			entropy,lyapunov = network_measures(stn,ensemble,N_steps)
		end
	end
			
	if return_stn
		return stn,retcode,entropy,lyapunov
	else
		return entropy,lyapunov
	end
		
end

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



