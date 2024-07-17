# StateTransitionNetworks.jl
Toolkit for dynamics on state-transition networks based on the paper Sándor et al. (2022).

## Installation
Download and install automatically to your local/global `julia` environment with:
```julia
] add https://github.com/rusandris/StateTransitionNetworks.jl.git
using StateTransitionNetworks
```

Or download a local copy, go to the directory and use the package in REPL:
```julia
] activate .
using StateTransitionNetworks
```

<!---

## Basic usage
### Discrete example: Henon map

Constructing state-transition network for the [Henon map](https://juliadynamics.github.io/DynamicalSystems.jl/dev/ds/predefined/#DynamicalSystemsBase.Systems.henon) (2D discrete dynamical system):
```julia
using DynamicalSystems
ds = PredefinedDynamicalSystems.henon()
traj, = trajectory(ds,10000;Ttr = 1000) #generate timeseries
```

Use `create_stn` to construct the graph object `stn` that corresponds to the STN:

```julia
stn,retcode = create_stn(traj,20); #using 20x20 grid
```
Internally, the phase space of the system is discretized using a grid of predefined size. The cells of this grid  will be the states (nodes) of our state-transition network. If the discretized trajectory contains a transition from state `i` to `j`, there will be an edge between two vertices `i -> j`.

The function also returns a `retcode` that informs the user about the state of the produced network which can take the following values:

* `:Success`  : graph is strongly connected
* `:NotConnected`  : graph isn't strongly connected
* `:Unusable`  : only when `make_ergodic` is enabled and the component with the most number of vertices is too small 


The network is a directed MetaGraph object which contains the vertices and edges of the STN which have metadata (attributes of predefined type) attached to them.

#### Vertex properties
* **Grid coordinates** : position (`:X` and `:y`) of the vertex point in the discretized phase space 
* **Probability **: probability of occurrence of the given state (vertex) 

```julia
stn[1] #returns the properties of vertex with label 1
Dict{Symbol, Union{Float64, Int64}} with 3 entries:
  :y    => 16
  :prob => 0.0354965
  :x    => 13
```

#### Edge properties 
* **Transition probabilities** (`:prob`) : conditional probability of the given transition 
* **Weights** (`:weight`) : probability of occurrence (non-conditional) of a given transition

```julia
stn[1,2] #returns properties of edge 1 -> 2
Dict{Symbol, Float64} with 2 entries:
  :weight => 0.0318
  :prob   => 0.895775
```

#### Stochastic (transition matrix) and weight matrix

Access the whole matrices `P[i,j]` and `Q[i,j]`

```julia
#sparse adjacency matrices 
P = get_transition_matrix(stn)
Q = get_weight_matrix(stn)
```

Use `network_measures` on the `stn` to calculate Sinai-Kolmogorov entropy and the Lyapunov network measure with numerical (random walk method) :
```julia
ensemble = 100 #number of random walks on the network
N_steps = 1e4 #number of steps taken in a random walk

entropy, lyapunov = network_measures(stn,ensemble, N_steps)
```

Use `network_measures` on the `P` stochastic transition matrix to calculate Sinai-Kolmogorov entropy and the Lyapunov network measure using the analytic formula :

```jul
entropy, lyapunov = network_measures(get_transition_matrix(stn))
```

Calculate the network measures for different parameters (dynamics) of the Henon map:

```julia
b = 0.3;
a_values = 1:0.001:1.4;
traj_length = 30000;
trans = 1000;
grid_size = 20;

entropy_measures = zeros(length(a_values))
lyapunov_measures = zeros(length(a_values))

for (i,a) in enumerate(a_values)
    system = PredefinedDynamicalSystems.henon([0.0, 0.0]; a=a, b=b)
    @show a
    timeseries, = trajectory(system, traj_length, [0, 0]; Ttr=trans)
    discrete_timeseries, vertex_names = timeseries_to_grid(timeseries, grid_size);
    stn, retcode = create_stn(timeseries, grid_size)
    entropy_measures[i], lyapunov_measures[i] = network_measures(get_transition_matrix(stn))
end
```

### Continuous example: Lorenz system
The workflow is similar for continuous systems. The difference is that one must turn the continuous system into map using [Poincaré surface of section](https://juliadynamics.github.io/DynamicalSystems.jl/dev/chaos/orbitdiagram/#Poincar%C3%A9-Surface-of-Section) (`PSOS`) with a given plane. This avoids the problem of high number of self-loops appearing in our STN.

Constructing state-transition network for the [Lorenz system](https://juliadynamics.github.io/DynamicalSystems.jl/dev/ds/predefined/#DynamicalSystemsBase.Systems.lorenz) (3D continuous dynamical system):
```julia
ds = PredefinedDynamicalSystems.lorenz()
plane = (1,15.0) #plane in 3D phase space with x = 15.0
n = 20 #using 20 x 20 grid
idxs = [2,3] #save y,z component of the PSOS
stn,retcode = create_stn(traj,n,plane,idxs)
```
## Time series analysis
There is a higher level function that accepts a time series and returns the corresponding `stn` and network measures.
```julia
stn_analysis(timeseries::Matrix;grid,plane,idxs,ensemble=100,N_steps=1000,make_ergodic=false, verbose=false,return_stn=false,use_analytic=false,use_stored_distribution=false)
```

-->


# References
Sándor, B.; Schneider, B.; Lázár, Z.I.; Ercsey-Ravasz M. : [A Novel Measure Inspired by Lyapunov Exponents for the Characterization of Dynamics in State-Transition Networks](https://www.mdpi.com/1099-4300/23/1/103)
