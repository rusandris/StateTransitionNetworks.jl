# StateTransitionNetworks.jl
Toolkit for dynamics on state-transition networks based on the paper Sándor et al. (2022).

## Basic usage

Using the package in REPL:
```
pkg> activate .
using StateTransitionNetworks
```
If one wants to return to the default environment to use other packages: `pkg> activate`

Constructing state-transition network for the Henon map (2D discrete dynamical system):
```julia
using DynamicalSystems
ds = Systems.henon()
traj = trajectory(ds,10000;Tr = 1500) #generate timeseries
```
This can be fed into `timeseries_to_grid` which discretizes the timeseries and returns the name of the vertices:
```
traj_grid, vertex_names = timeseries_to_grid(traj,20) # 20x20 grid
```
These vertices (nodes) will be the states of our state-transition network. If the discretized trajectory contains a transition from state 'i' to 'j', there will be an edge between two vertices `i -> j`.
Use `create_STN` to construct the graph object that corresponds to the STN:
```
stn_q, stn_p = create_STN(traj_grid,vertex_names)
```
The newtork is a SimpleWeightedDiGraph object. Two graphs are returned,one with occurence probability (`Q_ij`) and one with the transition probability
as weights (`P_ij`). 

Calculate entropy and lyapunov measures with `walk_statistics`:
```
ensemble = 100 #number of random walks on the network
N_steps = 1e4 #number of steps taken in a random walk

entropy, lyapunov = walk_statistics(ensemble, stn_p, N_steps)
```

### Generating figures: parameter dependence of the network measures 
Calculate the network measures for different parameters (dynamics) of the Henon map:

```
b = 0.3;
a_values = 1:0.001:1.4;
traj_length = 30000;
trans = 1000;
grid = 20;
ensemble = 100;
N_steps = 10000;

entropy = zeros(length(a_values))
lyapunov = zeros(length(a_values))

for (i,a) in enumerate(a_values)
    system = Systems.henon([0.0, 0.0]; a=a, b=b);
    timeseries = trajectory(system, traj_length, [0, 0]; Ttr=trans)
    discrete_timeseries, vertex_names = timeseries_to_grid(timeseries, grid);
    stn_q, stn_p = create_STN(discrete_timeseries, vertex_names)
    entropy[i], lyapunov[i] = walk_statistics(ensemble, stn_p, N_steps)
end
```
Plot results with:
```
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

# References
Sándor, B.; Schneider, B.; Lázár, Z.I.; Ercsey-Ravasz M. : [A Novel Measure Inspired by Lyapunov Exponents for the Characterization of Dynamics in State-Transition Networks](https://www.mdpi.com/1099-4300/23/1/103)
