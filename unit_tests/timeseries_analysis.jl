using StateTransitionNetworks
using DynamicalSystems
using Test


@testset "stn_analysis interface" begin

	ds = PredefinedDynamicalSystems.lorenz()	
	grid_size = 20
	plane1 = (1,15.0);
	plane2 = [1.0,0.0,0.0,0.0]
	idxs = [2,3]
	
	timeseries, = trajectory(ds,5000;Ttr = 1000,Δt = 1e-2)

	stn,retcode, nm1 = stn_analysis(timeseries; grid_size=grid_size,
		plane = plane1, 
		idxs = idxs, 
		ensemble = 100, 
		N_steps = 100, 
		make_ergodic = false, 
		verbose = false, 
		return_stn = true, 
		use_analytic = false, 
		use_stored_distribution = false,
		direction=-1)

	@test retcode == :Success
	
	nm2 = stn_analysis(timeseries; grid_size=grid_size,
		plane = plane2, 
		idxs = idxs,  
		make_ergodic = true, 
		verbose = false, 
		return_stn = false, 
		use_analytic = true, 
		use_stored_distribution = false)

	nm3 = stn_analysis(timeseries; grid_size=grid_size,
		plane = plane1, 
		idxs = idxs, 
		ensemble = 100, 
		N_steps = 100, 
		make_ergodic = true, 
		verbose = false, 
		return_stn = false, 
		use_analytic = true, 
		use_stored_distribution = true)

	@test nm1 != nm2 != nm3
end
