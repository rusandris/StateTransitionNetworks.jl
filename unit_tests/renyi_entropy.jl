using Revise
using StateTransitionNetworks
using DynamicalSystems
using Test

@testset "test truncated renyi_entropy at q=0 (topological entropy)" begin
    grid_size::Int = 20
    order::Int = 10 

    ds_logi = Systems.logistic() #r=4
    ds_henon = Systems.henon() #a=1.4, b=0.3
    orbit_logi,t = trajectory(ds_logi,10^6;Ttr=10^3)
    orbit_henon,t = trajectory(ds_henon,10^6;Ttr=10^3)
    
    sts_logi = zeros(Int128,length(orbit_logi))
    timeseries_to_grid!(sts_logi,orbit_logi, grid_size)
    sts_logi_higher_order = deepcopy(sts_logi)
    higher_order_symbolics!(sts_logi_higher_order,order)

    sts_henon = zeros(Int128,size(orbit_henon)[1])
    timeseries_to_grid!(sts_henon,orbit_henon, grid_size)
    sts_henon_higher_order = deepcopy(sts_henon)
    higher_order_symbolics!(sts_henon_higher_order,order)

    P,_,x =  calculate_transition_matrix(sts_logi)

    @testset "renyi_entropy" begin
        ht_logi = log(2) #topological entropy for logistic map (r=4)
        ht_henon = 0.4646 #topological entropy for henon map (a=1.4,b=0.3)
        P_logi_ho,_,_ = calculate_transition_matrix(@view sts_logi_higher_order[1:end-10])
        P_henon_ho,_,_ = calculate_transition_matrix(@view sts_henon_higher_order[1:end-10])
        K0_logi = renyi_entropy(P_logi_ho,0.0)
        K0_henon = renyi_entropy(P_henon_ho,0.0)
        @test isapprox(K0_logi,ht_logi;atol=1e-2)
        @test isapprox(K0_henon,ht_henon;atol=1e-2)

        #n is other than Inf
        renyi_entropy(P,1.2,1)

    end

    @testset "renyi_entropy_spectrums" begin    
        qs = [0.0:0.01:2;]
        renyi_entropy_spectrum(P,qs;x=x,verbose=false)
    end
end