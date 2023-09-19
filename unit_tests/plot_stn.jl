using StateTransitionNetworks
using DynamicalSystems
using Test

@testset "STN plotting" begin

	ds = PredefinedDynamicalSystems.henon()
	timeseries, = trajectory(ds,30000;Ttr=1000)
	
	stn,retcode = create_stn(timeseries,20;make_ergodic=false,verbose=false)
	plot_stn(stn;filename="plot_stn_test.pdf",nodesize=5, nodefillc="orange", linetype="curve", max_edgelinewidth=0.4, weight_exponent=1, nodelabels=true)
	
end

@testset "STN plotting: added timeseries" begin

	ds = PredefinedDynamicalSystems.lorenz()
	traj1, = trajectory(ds,1000;Ttr = 1000)
	traj2, = trajectory(ds,1000;Ttr = 1100)
	plane = (1,15.0)
	grid_size = 20
	
	stn12,retcode = add_timeseries([traj1,traj2],grid_size,plane;idxs=[2,3],verbose=true,make_ergodic=true,direction=+1) 
	
	plot_stn(stn12;filename="plot_stn_added_test.pdf",nodesize=5, nodefillc="orange", linetype="curve", max_edgelinewidth=0.4, weight_exponent=1, nodelabels=true)
	
end
