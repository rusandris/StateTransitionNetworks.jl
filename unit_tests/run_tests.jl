using Test

@testset verbose=true "StateTransitionNetworks" begin
	@testset verbose=true "gridding" include("timeseries_to_grid.jl")
	@testset verbose=true "transition_matrix" include("transition_matrix.jl")
	@testset verbose=true "network measures" include("network_measures.jl")
	@testset verbose=true "linsolve methods" include("iterative_linsolve.jl")
	@testset verbose=true "higher order states" include("higher_order_states.jl")
	@testset verbose=true "renyi entropy" include("renyi_entropy.jl")
	@testset verbose=true "SymbolStatistics structs" include("SymbolStatistics.jl")

		
	
	#old and deprecated (create_stn,stn plotting should work in the future, PRs needed)
	#@testset verbose=true "create_stn" include("create_stn.jl")
	#@testset verbose=true "stn plotting" include("plot_stn.jl")
	#@testset verbose=true "time series analysis" include("timeseries_analysis.jl")
	#@testset verbose=true "operations" include("operations.jl")
		
end;
