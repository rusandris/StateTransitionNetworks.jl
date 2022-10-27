"""
	stn_analysis(timeseries::AbstractDataset;grid,plane,ensemble,N_steps) ->  S, Λ
	
Calculates the Sinai-Kolmogorov Entropy and Lyapunov measure of a STN created from `timeseries`.
## Keyword arguments
* `grid` : size of grid used for discretization 
* `plane` : PSOS plane propagated to `poincaresos` from ChaosTools
* `ensemble` : number of individual random walks on the STN
* `N_steps` : number of steps in each random walk
"""
function stn_analysis(timeseries::AbstractDataset;grid,plane,ensemble=100,N_steps=10000)
	dim = size(timeseries)[2]
	psection = poincaresos(timeseries, plane); 
	psection = psection[:,filter(x -> x != plane[1],1:dim)]
	
	discrete_timeseries, vertex_names = timeseries_to_grid(psection,grid)
	stn = create_stn(discrete_timeseries,vertex_names)
	network_measures(stn,ensemble,N_steps)
end
