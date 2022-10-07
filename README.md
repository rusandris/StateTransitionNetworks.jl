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
traj = trajectory(ds,10000;Ttr = 1500) #generate timeseries
```
This can be fed into `timeseries_to_grid` which discretizes the timeseries and returns the name of the vertices:
```julia
traj_grid, vertex_names = timeseries_to_grid(traj,20) # 20x20 grid
```
These vertices (nodes) will be the states of our state-transition network. If the discretized trajectory contains a transition from state `i` to `j`, there will be an edge between two vertices `i -> j`.
Use `create_STN` to construct the graph object that corresponds to the STN:
```julia
stn_q, stn_p = create_STN(traj_grid,vertex_names)
```
The network is a SimpleWeightedDiGraph object. Two graphs are returned with a different set of weights: one with occurence probability weights (`Q_ij`) and one with the transition probability weights (`P_ij`). 

Calculate entropy and lyapunov measures with `walk_statistics`:
```julia
ensemble = 100 #number of random walks on the network
N_steps = 1e4 #number of steps taken in a random walk

entropy, lyapunov = walk_statistics(ensemble, stn_p, N_steps)
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
    timeseries = trajectory(system, traj_length, [0, 0]; Ttr=trans)
    discrete_timeseries, vertex_names = timeseries_to_grid(timeseries, grid_size);
    stn_q, stn_p = create_STN(discrete_timeseries, vertex_names)
    entropy_measures[i], lyapunov_measures[i] = walk_statistics(ensemble, stn_p, N_steps)
end
```
Plot results with:
```julia
using Plots, LaTeXStrings
plot(
   a_values, lyapunov,
   #ylims = (0,1.5),
   xticks = 1:0.05:1.4,
   linewidth = 4,
   linecolor = :black,
   xlabel = "a",
   ylabel = L"$\Lambda$",
   legend=false,
   xguidefontsize=18,
   yguidefontsize=18,
   tickfontsize=10
   )

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
stn_q, stn_p = create_STN(traj_grid,vertex_names)
```
We can also study the network measures for different parameters. With `PSOS`, the resulting network is not necesarrily strongly connected (not every vertex is reachable from every other vertex) in which case a random walk process would fail (we need ergodicity for the network measures). 

```julia
rho_values = 180:0.003:182;
T = 5000
grid_size = 20
ensemble = 100
N_steps = 10000
lyap_measures = zeros(length(rho_values))
entropy_measures = zeros(length(rho_values))
for (i,ρ) in enumerate(rho_values)
    ds = Systems.lorenz(ρ=ρ)
    @show ρ
    psection = poincaresos(ds, plane, T; Ttr=500, direction=+1,rootkw = (xrtol = 1e-8, atol = 1e-8))
    timeseries = psection[:,2:end]
    discrete_timeseries, vertex_names = timeseries_to_grid(timeseries, grid_size)
    _, stn_p = create_STN(discrete_timeseries, vertex_names)
    entropy_measures[i], lyap_measures[i] = walk_statistics(ensemble, stn_p, N_steps)
end
```

# References
Sándor, B.; Schneider, B.; Lázár, Z.I.; Ercsey-Ravasz M. : [A Novel Measure Inspired by Lyapunov Exponents for the Characterization of Dynamics in State-Transition Networks](https://www.mdpi.com/1099-4300/23/1/103)
