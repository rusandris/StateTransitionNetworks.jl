using StateTransitionNetworks
using DynamicalSystems
using Test

@testset "STN plotting" begin

	ds = PredefinedDynamicalSystems.henon()
	timeseries, = trajectory(ds,30000;Ttr=1000)
	
	stn,retcode = create_stn(timeseries,20;make_ergodic=false,verbose=false)
	plot_stn(stn;filename="plot_stn_test.pdf",nodesize=5, nodefillc="orange", linetype="curve", max_edgelinewidth=0.4, weight_exponent=1, nodelabels=true)
	
end
