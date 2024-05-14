using Revise
using StateTransitionNetworks
using DynamicalSystems
using Test

@testset "renyi_entropy" begin
    ds = Systems.henon()
    orbit,t = trajectory(ds,10^6;Ttr=10^3)
    grid_size::Int = 10
    sts = timeseries_to_grid(orbit, grid_size)
    P,_,x = calculate_transition_matrix(sts)
    
    renyi_entropy(P,1.2)
    renyi_entropy(P,0.0)

    #n is other than Inf
    renyi_entropy(P,1.2,1)

end

@testset "renyi_entropy_spectrums" begin
    ds = Systems.henon()
    orbit,t = trajectory(ds,10^6;Ttr=10^3)
    grid_size::Int = 10
    sts = timeseries_to_grid(orbit, grid_size)
    P,_,x = calculate_transition_matrix(sts)
    
    qs = [0.0:0.01:2;]
    renyi_entropy_spectrum(P,qs;x=x,verbose=false)

end