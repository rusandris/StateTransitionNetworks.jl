using StateTransitionNetworks
using DynamicalSystems
using Test


@testset "stn creation from stochastic matrix" begin
	P = [0 1. 1.;1. 0. 0.;0 0 0.]
	
	stn,retcode = create_stn(P)
	@test retcode == :NotConnected

	stn,retcode = create_stn(P;make_ergodic=true,verbose=true)
	@test retcode == :Success
	@test isnormalized(prob_matrix(stn))
	
end

@testset "stn creation from discrete timeseries" begin
	ds = PredefinedDynamicalSystems.henon()
	timeseries, = trajectory(ds,30000;Ttr=1000)
	incomplete_timeseries, = trajectory(ds,300;Ttr=0)
	
	stn,retcode = create_stn(timeseries,20;make_ergodic=false,verbose=false)
	@test retcode == :Success
	@test isnormalized(prob_matrix(stn))	
	
	inc_stn,inc_retcode = create_stn(incomplete_timeseries,20;make_ergodic=false,verbose=true)
	@test inc_retcode == :NotConnected
	
	inc_stn,inc_retcode = create_stn(incomplete_timeseries,20;make_ergodic=true,verbose=true)
	@test inc_retcode == :Success
	@test isnormalized(prob_matrix(inc_stn))
end


@testset "stn creation from continuous timeseries" begin
	ds = PredefinedDynamicalSystems.lorenz()	
	grid_size = 20
	plane1 = (1,15.0);
	plane2 = [1.0,0.0,0.0,0.0]
	idxs = [2,3]
	
	
	timeseries, = trajectory(ds,5000;Ttr = 1000,Δt = 1e-2)
	incomplete_timeseries, = trajectory(ds,50;Ttr = 0,Δt = 1e-2)

    stn, retcode = create_stn(timeseries, grid_size,plane1,idxs; make_ergodic=false,verbose=false);
    @test retcode == :Success
    @test isnormalized(prob_matrix(stn))
    
    stn, retcode = create_stn(timeseries, grid_size,plane2,idxs; make_ergodic=false,verbose=false);
    @test retcode == :Success
    @test isnormalized(prob_matrix(stn))
    
    inc_stn, inc_retcode = create_stn(incomplete_timeseries, grid_size,plane1,idxs; make_ergodic=false,verbose=false);
    @test inc_retcode == :NotConnected
    
    inc_stn, inc_retcode = create_stn(incomplete_timeseries, grid_size,plane1,idxs; make_ergodic=true,verbose=false);
    @test inc_retcode == :Unusable
	
	
end

@testset "miscallenous API" begin
	ds = PredefinedDynamicalSystems.lorenz()	
	grid_size = 20
	plane1 = (1,15.0);
	idxs = [2,3]	
	timeseries, = trajectory(ds,5000;Ttr = 1000,Δt = 1e-2)
    stn, retcode = create_stn(timeseries, grid_size,plane1,idxs; make_ergodic=false,verbose=false);

	
	state_distrib,state_pos = state_distribution(stn)
	@test sum(state_distrib) ≈ 1.0
	@test isnormalized(prob_matrix(stn))
	@test sum(weight_matrix(stn))  ≈ 1.0
	renormalize!(stn)
	@test isnormalized(prob_matrix(stn))
	@test sum(weight_matrix(stn))  ≈ 1.0
	

end



