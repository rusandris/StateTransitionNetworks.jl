using Revise
using StateTransitionNetworks
using SparseArrays
using DynamicalSystems
using Random
using Test


@testset "Stochasticity and ergodicity" begin
    P = sparse([1.0 1.0 0; 1.0 1.0 0; 0 0 1.0])
    @test !(is_stochastic(P))
    @test !(is_strongly_connected(P))

    Pn = calculate_transition_matrix(P) #normalize
    @test is_stochastic(Pn)

    Pc = sparse([1.0 1.0; 1.0 1.0 ]) #make it connected
    @test is_strongly_connected(Pc)
end

@testset "2D system tests" begin

    ds = Systems.henon()
    orbit,t = trajectory(ds,10^6;Ttr=10^3)
    grid_size::Int = 10
    sts = timeseries_to_grid(orbit, grid_size)
    symbs = unique(sts)
    nr_symbs = length(symbs)

    @testset "calculate_transition_matrix (oop) and (ip)" begin
       

        P1,Q,x = calculate_transition_matrix(sts)

        @test is_stochastic(P1)
        @test is_strongly_connected(P1)
        @test sum(Q) ≈ 1     
        @test sum(x) ≈ 1
        
        symbol_dictionary = Dict(shuffle(symbs) .=> 1:nr_symbs)
        P2,_,_ = calculate_transition_matrix(sts;symbol_dictionary=symbol_dictionary,verbose=true)
        @test is_stochastic(P2)
        @test is_strongly_connected(P2)

        #test if we still get the same measures even if the matrix is shuffled
        @test all(network_measures(P1) .≈ network_measures(P2))

        sts_copy = deepcopy(sts)
        P,Q,x = calculate_transition_matrix!(sts_copy)

        #test if sts_copy was remapped in-place
        @test !(all(sts .== sts_copy))
        @test all(unique(sts_copy) .== 1:nr_symbs)

        @test is_stochastic(P)
        @test is_strongly_connected(P)
        @test sum(Q) ≈ 1     
        @test sum(x) ≈ 1
    

        #test if we still get the same measures even if the matrix is shuffled
        @test all(network_measures(P1) .≈ network_measures(P2))
    end

end

@testset "1D system tests" begin

    ds = Systems.logistic()
    orbit,t = trajectory(ds,10^6;Ttr=10^3)
    grid_size::Int = 32
    sts = timeseries_to_grid(orbit, grid_size)
    symbs = unique(sts)
    nr_symbs = length(symbs)

    @testset "calculate_transition_matrix_noremap" begin

        #for the logistic map (r=4)
        #it can be assumed that all symbols appear
        #and they go from 1:32
        
        P,Q,x = calculate_transition_matrix_no_remap(sts)

        #test if sts_copy was remapped in-place

        @test is_stochastic(P)
        @test is_strongly_connected(P)
        @test sum(Q) ≈ 1     
        @test sum(x) ≈ 1
    end

    @testset "calculate_transition_matrix (oop) and (ip)" begin
       

        P1,Q,x = calculate_transition_matrix(sts)

        @test is_stochastic(P1)
        @test is_strongly_connected(P1)
        @test sum(Q) ≈ 1     
        @test sum(x) ≈ 1
        
        symbol_dictionary = Dict(shuffle(symbs) .=> 1:nr_symbs)
        P2,_,_ = calculate_transition_matrix(sts;symbol_dictionary=symbol_dictionary,verbose=true)

        #test if we still get the same measures even if the matrix is shuffled
        @test all(network_measures(P1) .≈ network_measures(P2))

        sts_copy = deepcopy(sts)
        P,Q,x = calculate_transition_matrix!(sts_copy)

        #test if sts_copy was remapped in-place
        @test !(all(sts .== sts_copy))
        @test all(unique(sts_copy) .== 1:nr_symbs)

        @test is_stochastic(P)
        @test is_strongly_connected(P)
        @test sum(Q) ≈ 1     
        @test sum(x) ≈ 1
    

        #test if we still get the same measures even if the matrix is shuffled
        @test all(network_measures(P1) .≈ network_measures(P2))
    end

end


