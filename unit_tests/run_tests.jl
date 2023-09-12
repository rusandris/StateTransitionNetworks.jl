using StateTransitionNetworks
using DynamicalSystems
using Test

@testset verbose=true "StateTransitionNetworks" begin
	@testset verbose=true "create_stn" include("create_stn.jl") 
	@testset verbose=true "operations" include("operations.jl")
	@testset verbose=true "network measures" include("network_measures.jl")	
	@testset verbose=true "stn plotting" include("plot_stn.jl")	
	@testset verbose=true "time series analysis" include("timeseries_analysis.jl")	
end
;
