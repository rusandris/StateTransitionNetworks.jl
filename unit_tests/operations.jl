using StateTransitionNetworks
using DynamicalSystems
using Test


@testset "stn addition test: addition" begin
	stn_added_test,retcode = add_timeseries([[(1,2),(2,3),(3,4),(1,2),(3,4)],[(3,4),(10,4),(1,2)]],10)
	@test retcode == :Success
	@test length(prob_matrix(stn_added_test).nzval) == 6
end


ds = PredefinedDynamicalSystems.henon()
traj1, = trajectory(ds,10000;Ttr = 1000)
traj2, = trajectory(ds,10000;Ttr = 11000)


@testset "stn addition test: commutativity" begin
	symts_list12 = timeseries_to_common_grid([traj1,traj2],20)
	symts_list21 = timeseries_to_common_grid([traj2,traj1],20)
	
	stn_added12,retcode = add_timeseries(symts_list12,20)
	stn_added21,retcode = add_timeseries(symts_list21,20)
	@test are_equal(stn_added12,stn_added21)
	
end
