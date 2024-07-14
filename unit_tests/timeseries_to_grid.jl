using StateTransitionNetworks
using DynamicalSystems
using Test



#---------------------------------2D grid tests----------------------------------

@testset "grid partitioning of 2D random points" begin
    grid_size = 5
    points = rand(10,2)
    points_bad = deepcopy(points)
    points_bad[end,:] .= [100.0,100.0] #bad, large value
    grid_edges=[0.0,0.0,1.0,1.0]

    #-----------------------------------------out-of-place---------------------------------------

    sts0 = timeseries_to_grid(points,grid_size) #tries to include everything
    sts01 = timeseries_to_grid(points_bad,grid_size) #tries to include everything

    @test_throws DomainError timeseries_to_grid(points_bad,grid_size,grid_edges=grid_edges,outside_grid=:error) #errors if outside
    sts1 = timeseries_to_grid(points,grid_size,grid_edges=grid_edges) #without outsiders
    sts2 = timeseries_to_grid(points_bad,grid_size,grid_edges=grid_edges,outside_grid=:skip) #skips if outside
    sts3 = timeseries_to_grid(points_bad,grid_size,grid_edges=grid_edges,outside_grid=:include) #includes artificially if outside
    
    #test for the lengths of the symbolic ts
    @test length(sts1) == 10
    @test length(sts2) == 9
    @test length(sts3) == 10

    #test for the values of the symbolic ts
    @test all(x -> x in 1:grid_size^2,sts1)
    @test all(x -> x in 1:grid_size^2,sts2)
    @test all(x -> x in 1:grid_size^2,sts3)

    #test if symbols are the same except the bad one
    @test all(sts1[1:9] .== sts3[1:9])
    @test all(sts2 .== sts1[1:9])
    @test all(sts2 .== sts3[1:9])


    #-----------------------------------------in-place---------------------------------------

    sts = zeros(Int,10)

    timeseries_to_grid!(sts,points_bad,grid_size) #normal
    @test_throws ArgumentError timeseries_to_grid!(sts,points_bad,grid_size,grid_edges=grid_edges,outside_grid=:skip) #not an option
    @test_throws DomainError timeseries_to_grid!(sts,points_bad,grid_size,grid_edges=grid_edges,outside_grid=:error) #skips if outside
    timeseries_to_grid!(sts,points_bad,grid_size,grid_edges=grid_edges,outside_grid=:include) #includes artificially if outside

end

@testset "grid partitioning on 2D orbit" begin
    grid_size=20
    grid_edges=[0.0,0.0,0.5,0.5]

    ds = Systems.henon()
    orbit,t = trajectory(ds,10^7,[0.4,0.4];Ttr=1000)
    sts1 = timeseries_to_grid(orbit,grid_size) #without outsiders
    sts2 = timeseries_to_grid(orbit,grid_size,grid_edges=grid_edges,outside_grid=:skip) #with skipping outsiders

end

#-----------------------------------------------1D grid tests-----------------------------------------------

@testset "grid partitioning of 1D random points" begin
    grid_size = 5
    points = rand(10)
    points_bad = deepcopy(points)
    points_bad[end] = 100.0 #bad, large value
    grid_edges=[0.0,1.0]

    #-----------------------------------------out-of-place---------------------------------------

    sts0 = timeseries_to_grid(points,grid_size) #tries to include everything
    sts01 = timeseries_to_grid(points_bad,grid_size) #tries to include everything

    @test_throws DomainError timeseries_to_grid(points_bad,grid_size,grid_edges=grid_edges,outside_grid=:error) #errors if outside
    sts1 = timeseries_to_grid(points,grid_size,grid_edges=grid_edges) #without outsiders
    sts2 = timeseries_to_grid(points_bad,grid_size,grid_edges=grid_edges,outside_grid=:skip) #skips if outside
    sts3 = timeseries_to_grid(points_bad,grid_size,grid_edges=grid_edges,outside_grid=:include) #includes artificially if outside
    
    #test for the lengths of the symbolic ts
    @test length(sts1) == 10
    @test length(sts2) == 9
    @test length(sts3) == 10

    #test for the values of the symbolic ts
    @test all(x -> x in 1:grid_size^2,sts1)
    @test all(x -> x in 1:grid_size^2,sts2)
    @test all(x -> x in 1:grid_size^2,sts3)

    #test if symbols are the same except the bad one
    @test all(sts1[1:9] .== sts3[1:9])
    @test all(sts2 .== sts1[1:9])
    @test all(sts2 .== sts3[1:9])


    #-----------------------------------------in-place---------------------------------------

    sts = zeros(Int,10)

    timeseries_to_grid!(sts,points_bad,grid_size) #normal
    @test_throws ArgumentError timeseries_to_grid!(sts,points_bad,grid_size,grid_edges=grid_edges,outside_grid=:skip) #not an option
    @test_throws DomainError timeseries_to_grid!(sts,points_bad,grid_size,grid_edges=grid_edges,outside_grid=:error) #skips if outside
    timeseries_to_grid!(sts,points_bad,grid_size,grid_edges=grid_edges,outside_grid=:include) #includes artificially if outside

end

@testset "grid partitioning on 1D orbit" begin
    grid_size=20
    grid_edges_bad = [0.0,0.0,0.5,0.5]
    grid_edges = [0.0,0.5]

    ds = Systems.logistic()
    orbit,t = trajectory(ds,10^7,[0.4];Ttr=1000)
    sts1 = timeseries_to_grid(orbit,grid_size) #without outsiders
    @test_throws ArgumentError timeseries_to_grid(orbit,grid_size,grid_edges=grid_edges_bad) #with skipping outsiders
    sts2 = timeseries_to_grid(orbit,grid_size,grid_edges=grid_edges,outside_grid=:skip) #with skipping outsiders
    @test length(sts1) != length(sts2)
end


#-----------------------------------------------general 1D grid tests-----------------------------------------------
#generating partition of the tent map

@testset "generating partition of the tent map" begin

    #tent map definition
    function f(u, p, n)
        r = p[1]
        x = u[1]
        if x<=r
        return SVector(x/r)
        else
        return SVector((1-x)/(1-r))
        end
    end
    
    r = 0.6
    ds = DeterministicIteratedMap(f, [0.4], [r])
    orbit,t = trajectory(ds,10^4;Ttr=1000)

    tent_generating_partition(x::Float64,r::Float64) = x <= r ? 1 : 2

    sts = timeseries_to_grid(orbit,x -> tent_generating_partition(x,r))
    @test length(unique(sts)) == 2

end