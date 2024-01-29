using StateTransitionNetworks
using DynamicalSystems
using Test


@testset "network measures" begin

	as = [1.2,1.2265,1.27,1.24]
	
	ensemble = 10000
	nr_steps = 10000

	for (i,a) in enumerate(as)
		ds = PredefinedDynamicalSystems.henon(a=a)
		timeseries, = trajectory(ds,30000;Ttr=1000)
		stn,retcode = create_stn(timeseries,20;make_ergodic=false,verbose=false)
		
		#numerical (random walk) method
		num_nm = network_measures(stn,ensemble,nr_steps)
		
		@show a
		sts = timeseries_to_grid(timeseries,20)
		P = calculate_transition_matrix(sts)
		#analytic formula
		analytic_nm = network_measures(P)
		#lyapunov_measure(P;maxiter=1000)
		
		#analytic formula - presaved stationary distribution
		presaved_distrib = state_distribution(stn)[1]
		analytic_presaved_nm = network_measures(prob_matrix(stn);x=presaved_distrib)		
		
		if i == 3 #periodic
			@test all(x -> x ≈ 0,num_nm)
			@test all(x -> x ≈ 0,analytic_nm)

		else #chaotic
			@test all(x -> x > 0, num_nm)
			@test all(x -> x > 0,analytic_nm)

		end
		
		@test all(isapprox.(analytic_nm,num_nm,atol=1e-1))
		@test all(isapprox.(analytic_nm,analytic_presaved_nm,atol=1e-4))
		
	end
	

end
