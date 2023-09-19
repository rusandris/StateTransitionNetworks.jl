using StateTransitionNetworks
using DynamicalSystems
using Test


@testset "stn addition test: addition" begin
	stn_added_test,retcode = add_timeseries([[(1,2),(2,3),(3,4),(1,2),(3,4)],[(3,4),(10,4),(1,2)]],10)
	@test retcode == :Success
	@test length(prob_matrix(stn_added_test).nzval) == 6
end



@testset "stn addition test: commutativity" begin
	ds = PredefinedDynamicalSystems.henon()
	traj1, = trajectory(ds,10000;Ttr = 1000)
	traj2, = trajectory(ds,10000;Ttr = 11000)

	symts_list12 = timeseries_to_common_grid([traj1,traj2],20)
	symts_list21 = timeseries_to_common_grid([traj2,traj1],20)
	
	stn_added12,retcode = add_timeseries(symts_list12,20)
	stn_added21,retcode = add_timeseries(symts_list21,20)
	@test are_equal(stn_added12,stn_added21)
	
end

@testset "stn addition test: continuous timeseries" begin

	ds = PredefinedDynamicalSystems.lorenz()
	traj1, = trajectory(ds,1000;Ttr = 1000)
	traj2, = trajectory(ds,1000;Ttr = 1100)
	plane = (1,15.0)
	grid_size = 20
	
	stn12,retcode = add_timeseries([traj1,traj2],grid_size,plane;idxs=[2,3],verbose=true,make_ergodic=true,direction=+1) 
	stn21,retcode = add_timeseries([traj2,traj1],grid_size,plane;idxs=[2,3],verbose=true,make_ergodic=true,direction=+1) 
	
	state_distrib12,state_pos12 = state_distribution(stn12) 
	state_distrib21,state_pos21 = state_distribution(stn21)
	
	@test are_equal(stn12,stn21)
	#@test state_distrib12 == state_distrib21
	#@test state_pos12 == state_pos21

end
