# StateTransitionNetworks.jl
Toolkit for dynamics on state-transition networks based on the paper Sándor et al. (2022).


## Installation
Download a local copy, go to the directory and use the package in REPL:
```julia
pkg> activate .
using StateTransitionNetworks
```
Or download and install automatically to your local/global `julia` environment with:
```julia
pkg> add https://github.com/rusandris/StateTransitionNetworks.jl.git
using StateTransitionNetworks
```


## Basic usage
### Discrete example: Henon map

Constructing state-transition network for the [Henon map](https://juliadynamics.github.io/DynamicalSystems.jl/dev/ds/predefined/#DynamicalSystemsBase.Systems.henon) (2D discrete dynamical system):
```julia
using DynamicalSystems
ds = Systems.henon()
traj = trajectory(ds,10000;Ttr = 1000) #generate timeseries
```
This can be fed into `timeseries_to_grid` which discretizes the timeseries and returns the name of the vertices:
```julia
traj_grid, vertex_names = timeseries_to_grid(traj,20) # 20x20 grid
```
These vertices (nodes) will be the states of our state-transition network. If the discretized trajectory contains a transition from state `i` to `j`, there will be an edge between two vertices `i -> j`.
Use `create_stn` to construct the graph object `stn` that corresponds to the STN:
```julia
stn,retcode = create_stn(traj_grid,vertex_names); #output might be long
```
The function also returns a `retcode` which can take the following values:
* `:Success` : every vertex has at least one incoming and outgoing edge
* `:NoIncoming` : the graph contains vertex(es) that have no incoming edges
* `:NoOutgoing` : the graph contains vertex(es) that have no outgoing edges
* `:DeadEnd` : self-loop vertex (fixed point)


The network is a directed MetaGraph object which contains the vertices and edges of the STN which have metadata (attributes of predefined type) attached to them.

```julia

 MetaGraph(DiGraph(),Label = Int64, 
		VertexData = Dict{Symbol, Int64}, 
		EdgeData = Dict{Symbol,Float64}, 
		default_weight = 0.0)
```

#### Vertex properties
* **Grid coordinates** (position of the vertex point in the discretized phase space) 
#### Edge properties 
* **Transition probabilities** and **weights**

Accesing vertex and edge properties:
```julia
stn[1] #returns properties of vertex with label 1
Dict{Symbol, Int64} with 2 entries:
  :y => 15
  :x => 12

stn[1,2] #returns properties of edge 1 -> 2
Dict{Symbol, Float64} with 2 entries:
  :weight => 0.023
  :prob   => 0.741935

```


Access the whole matrices `P[i,j]` and `Q[i,j]`
```julia
#sparse adjacency matrices 
P = prob_matrix(stn)
Q = weight_matrix(stn)
```

Calculate entropy and lyapunov measures with `network_measures`:
```julia
ensemble = 100 #number of random walks on the network
N_steps = 1e4 #number of steps taken in a random walk

entropy, lyapunov = network_measures(stn,ensemble, N_steps)
```



Calculate the network measures for different parameters (dynamics) of the Henon map:

```julia
b = 0.3;
a_values = 1:0.001:1.4;
traj_length = 30000;
trans = 1000;
grid_size = 20;
ensemble = 100;
N_steps = 10000;

entropy_measures = zeros(length(a_values))
lyapunov_measures = zeros(length(a_values))

for (i,a) in enumerate(a_values)
    system = Systems.henon([0.0, 0.0]; a=a, b=b)
    @show a
    timeseries,  = trajectory(system, traj_length, [0, 0]; Ttr=trans)
    discrete_timeseries, vertex_names = timeseries_to_grid(timeseries, grid_size);
    stn,retcode = create_stn(discrete_timeseries, vertex_names)
    entropy_measures[i], lyapunov_measures[i] = network_measures(stn,ensemble, N_steps)
end
```

### Continuous example: Lorenz system
The workflow is similar for continuous systems. The difference is that one must turn the continuous system into map using [Poincaré surface of section](https://juliadynamics.github.io/DynamicalSystems.jl/dev/chaos/orbitdiagram/#Poincar%C3%A9-Surface-of-Section) (`PSOS`) with a given plane. This avoids the problem of high number of self-loops appearing in our STN.

Constructing state-transition network for the [Lorenz system](https://juliadynamics.github.io/DynamicalSystems.jl/dev/ds/predefined/#DynamicalSystemsBase.Systems.lorenz) (3D continuous dynamical system):
```julia
ds = Systems.lorenz()
plane = (1,15.0) #plane in 3D phase space with x = 15.0
```
Calculate the `PSOS` with `poincaresos`

```julia
T = 500
psection = poincaresos(ds, plane, T; Ttr=300, direction=+1); 
traj = psection[:,2:end] #select y,z variables only
```

> **_NOTE:_** The default tolerances for the `poincaresos` are `rootkw = (xrtol = 1e-6, atol = 1e-6)`. These might not be low enough for all purposes (for ex. periodic attractors).

This can be treated as a 2D map:
```julia
traj_grid, vertex_names = timeseries_to_grid(traj,20) # 20x20 grid
stn,retcode = create_stn(traj_grid,vertex_names)
```
We can also study the network measures for different parameters. 
```julia
rho_values = 180:0.003:182;
T = 5000
grid_size = 20
ensemble = 100
N_steps = 10000
plane = (1,15.0) #plane in 3D phase space with x = 15.0
lyap_measures = zeros(length(rho_values))
entropy_measures = zeros(length(rho_values))
for (i,ρ) in enumerate(rho_values)
    ds = Systems.lorenz(ρ=ρ)
    @show ρ
    psection = DynamicalSystemsBase.poincaresos(ds, plane, T; Ttr=500, idxs = [2,3],direction=+1,rootkw = (xrtol = 1e-8, atol = 1e-8))
    discrete_timeseries, vertex_names = timeseries_to_grid(psection, grid_size)
    stn,retcode = create_stn(discrete_timeseries, vertex_names)
    entropy_measures[i], lyap_measures[i] = network_measures(stn,ensemble, N_steps)
end
```
## Network measures from data
There is a function that accepts a timeseries and returns the corresponding `stn` and network measures.
```julia
stn_analysis(timeseries::Matrix;grid,plane,idxs,ensemble=100,N_steps=1000,return_stn=false)
```

#### Keyword arguments
* `grid` : size of grid used for discretization 
* `plane` : PSOS plane propagated to `poincaresos` from [ChaosTools](https://github.com/JuliaDynamics/ChaosTools.jl)
* `idxs` : choose which variables to save
* `ensemble` : number of individual random walks on the STN
* `N_steps` : number of steps in each random walk
* `return_stn` : returns only the `stn` and the `retcode`






# References
Sándor, B.; Schneider, B.; Lázár, Z.I.; Ercsey-Ravasz M. : [A Novel Measure Inspired by Lyapunov Exponents for the Characterization of Dynamics in State-Transition Networks](https://www.mdpi.com/1099-4300/23/1/103)
