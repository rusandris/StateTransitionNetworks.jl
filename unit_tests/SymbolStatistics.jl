using Revise
using StateTransitionNetworks
using DynamicalSystems
using Random
using Test



@testset " SymbolStatistics structs: " begin 

    ds = Systems.henon()
    orbit,t = trajectory(ds,10^6;Ttr=10^3)
    grid_size::Int = 10
    sts = timeseries_to_grid(orbit, grid_size)

    @testset "Henon map " begin
        trs = SymbolStatistics(Int64, grid_size*grid_size)
        add_transitions!(trs, sts)
        P_trs = calculate_transition_matrix(trs)
        @test is_stochastic(P_trs)
        x_trs = calculate_symbol_probabilities(trs)
        @test sum(x_trs) .≈ 1
        P,Q,x = calculate_transition_matrix(sts)
        @test all(network_measures(P_trs) .≈ network_measures(P))             
    end
    
    @testset "Test transition addition" begin

        @testset "SymbolStatistics" begin
            #test adding single transitions
            trs = SymbolStatistics(Int64,3)
            add_transition!(trs, (11,22))
            @test trs.nr_unique_symbols == 2
            @test trs.nr_transitions == 1
            @test trs.last_transition == (11,22)

            add_transition!(trs, (22,33))
            @test trs.nr_unique_symbols == 3
            @test trs.nr_transitions == 2
            @test sum(trs.x) == 3

            #test if adding single transitions and adding transitions
            #from a symbolic trajectory is the same
            trs_in_one_go = SymbolStatistics(Int64,3)
            add_transitions!(trs_in_one_go,[11,22,33])
            @test all([getfield(trs,field) == getfield(trs,field) for field in fieldnames(typeof(trs))]) 
        end

        #=
        @testset "more general (string) symbols" begin
            trs = SymbolStatistics(String,10)
            add_transitions!(trs, ["a","b","c","b","a","c"])
            @test trs.nr_unique_symbols == 3
            @test trs.nr_transitions == 5     
            P_trs = calculate_transition_matrix(trs)
            @test is_stochastic(P_trs)
            x_trs = calculate_symbol_probabilities(trs)

            P,_,x = calculate_transition_matrix([1,2,3,2,1,3])
            @test all(P_trs .== P)
            @test all(network_measures(P_trs) .≈ network_measures(P))
        end
        =#
    end

    @testset "SymbolStatisticsNoRemap" begin
        #test adding single transitions
        trs = SymbolStatisticsNoRemap(Int64,3)
        add_transition!(trs, (1,2))
        @test trs.nr_unique_symbols == 2
        @test trs.nr_transitions == 1
        @test trs.last_transition == (1,2)

        add_transition!(trs, (2,3))
        @test trs.nr_unique_symbols == 3
        @test trs.nr_transitions == 2
        @test sum(trs.x) == 3

        #test if adding single transitions and adding transitions
        #from a symbolic trajectory is the same
        trs_in_one_go = SymbolStatisticsNoRemap(Int64,3)
        add_transitions!(trs_in_one_go,[1,2,3])
        @test all([getfield(trs,field) == getfield(trs,field) for field in fieldnames(typeof(trs))]) 
    end


end


