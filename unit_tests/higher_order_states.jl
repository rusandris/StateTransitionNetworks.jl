using Revise
using StateTransitionNetworks
using DynamicalSystems
using Random
using Test

@testset "higher order states" begin
    grid_size::Int = 100
    order_low::Int = 3
    order_high::Int = 10 

    ds_logi = Systems.logistic() #r=4
    orbit_logi,t = trajectory(ds_logi,10^6;Ttr=10^3)

    sts_logi = zeros(Int128,length(orbit_logi))
    timeseries_to_grid!(sts_logi,orbit_logi, grid_size)
    sts_logi_higher_order = deepcopy(sts_logi)
    higher_order_symbolics!(sts_logi_higher_order,order_high)

    sts_logi = timeseries_to_grid(orbit_logi, grid_size)
    #test if overflow happens
    @test_logs (:warn,) higher_order_symbolics!(sts_logi,order_high)
        


end