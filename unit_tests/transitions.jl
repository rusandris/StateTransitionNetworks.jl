using Revise
using StateTransitionNetworks
using DynamicalSystems
using Random
using Test

@testset " test transition adding: " begin 
    ds = Systems.henon()
    orbit,t = trajectory(ds,10^6;Ttr=10^3)
    grid_size::Int = 10
    sts = timeseries_to_grid(orbit, grid_size)

    @testset "IntegerTransitions " begin
        trs = IntegerTransitions(Int64, grid_size*grid_size)
        add_transitions!(trs, sts)
        P_trs = calculate_transition_matrix(trs)
        @test is_stochastic(P_trs)
        x_trs = calculate_symbol_probabilities(trs)
        @test sum(x_trs) .≈ 1
        P,Q,x = calculate_transition_matrix(sts)
        @test all(network_measures(P_trs) .≈ network_measures(P))             
    end

    @testset "GeneralTransitions " begin
        trs = GeneralTransitions(String,10)
        add_transitions!(trs, ["a","b","c","b","a","c"])
        @test trs.nr_used_symbols == 3
        @test trs.nr_transitions == 5     
        P_trs = calculate_transition_matrix(trs)
        @test is_stochastic(P_trs)
        x_trs = calculate_symbol_probabilities(trs)

        P,_,x = calculate_transition_matrix([1,2,3,2,1,3])
        @test all(P_trs .== P)
        @test all(network_measures(P_trs) .≈ network_measures(P))

    end

end