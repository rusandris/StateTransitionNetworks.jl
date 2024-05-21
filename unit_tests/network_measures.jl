using StateTransitionNetworks
using DynamicalSystems
using Test


@testset "test analytic formula and random walk sim" begin

	as = [1.2,1.2265,1.24,1.27]
	
	ensemble = 1000
	nr_steps = 1000

	for (i,a) in enumerate(as)
		ds = Systems.henon(a=a)
		timeseries, = trajectory(ds,30000;Ttr=1000)
		sts = timeseries_to_grid(timeseries,20)
		P,Q,x = calculate_transition_matrix(sts)
		
		#numerical (random walk) method
		num_nm = network_measures(P,ensemble,nr_steps)
		
		@show a
		
		#analytic formula
		analytic_nm = network_measures(P)
		@show analytic_nm
		#lyapunov_measure(P;maxiter=1000)
		
		analytic_presaved_nm = network_measures(P;x=x)		
		
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

@testset "test measure calculation options and methods" begin

	ds = Systems.logistic()
	orbit,t = trajectory(ds,10^6;Ttr=10^3)
	grid_size::Int = 10
	sts = timeseries_to_grid(orbit, grid_size)
	symbs = unique(sts)
	nr_symbs = length(symbs)

	P,Q,x = calculate_transition_matrix(sts)

	@testset "test Λ linalg solvers" begin


		ϵ = 1e-4
		maxiter = 10000
		S1,Λ1 = network_measures(P;ϵ=ϵ,maxiter=maxiter,alg=iterative_linsolve) #own iterative method
		S2,Λ2 = network_measures(P;ϵ=ϵ,maxiter=maxiter,alg=linsolve) #KrylovKit method
		S3,Λ3 = network_measures(P;ϵ=ϵ,maxiter=maxiter,alg=hybrid_solve) #combination of the above

		#@test S1 ≈ S2 ≈ S3
		@test isapprox(S1, S2; atol=0.01)
		@test isapprox(S2, S3; atol=0.01) 
		
	end

	@testset "bit_number_measures" begin
		C1,C2 = bit_number_measures(x)
		x_calc = stationary_distribution(P)
		C1_calc,C2_calc = bit_number_measures(x_calc)
		@test isapprox(C1,C1_calc;atol=1e-3)
		@test isapprox(C2,C2_calc;atol=1e-3)
	end

end