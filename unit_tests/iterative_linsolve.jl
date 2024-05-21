using StateTransitionNetworks
using DynamicalSystems
using Test
using LinearAlgebra


@testset "test linsolve algs" begin
    ds = Systems.logistic()
	orbit,t = trajectory(ds,10^6;Ttr=10^3)
	grid_size_low::Int = 10
    grid_size_high::Int = 2000
	sts_low = timeseries_to_grid(orbit, grid_size_low)
    sts_high = timeseries_to_grid(orbit, grid_size_high)


	P_small,Q_small,x_small = calculate_transition_matrix(sts_low)
    P_large,Q_large,x_large = calculate_transition_matrix(sts_high)

    X_small = PseudoDenseMatrix(x_small)
    L_small = -sparse_log(P_small)
    L_small = P_small .* L_small
    v_small = ones(length(x_small))

    X_large = PseudoDenseMatrix(x_large)
    L_large = -sparse_log(P_large)
    L_large = P_large .* L_large
    v_large = ones(length(x_large))

    ϵ = 1e-4
    maxiter_high = 10000
    maxiter_low = 2

    z1_large,conv_info = iterative_linsolve(P_large,X_large,L_large*v_large;ϵ = ϵ,maxiter=maxiter_high)
    z1_small,conv_info = iterative_linsolve(P_small,X_small,L_small*v_small;ϵ = ϵ,maxiter=maxiter_high) 
    z2_large,conv_info = linsolve(I - P_large + X_large,L_large*v_large;tol = ϵ,maxiter=maxiter_high)
    z2_small,conv_info = linsolve(I - P_small + X_small,L_small*v_small;tol = ϵ,maxiter=maxiter_high)

    #test if they give the same results
    @test all(isapprox.(z1_large,z2_large;atol=1e-3))
    @test all(isapprox.(z1_small,z2_small;atol=1e-3) )

    #test alg switching 
    #large matrix doesn't converge with iterative_linsolve -> switch to KrylovKit.linsolve 
    @test_logs (:warn,) hybrid_solve(P_large,X_large,L_large*v_large;ϵ = ϵ,maxiter=maxiter_low) 



end