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

## Basic usage
Given a time series, like
```julia
x = [0.1, 0.2, 0.7, 0.3, 0.1, 0.3, 0.2, 0.1, 0.3, 0.4, 0.1, 0.0] #example time series
```

#### Create symbols  (also called outcomes) from the data points of the time series 
##### For example using [ComplexityMeasures.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/dev/)
* by applying binning
```julia
using ComplexityMeasures
binning = ValueBinning(RectangularBinning(4)) #4 bins (symbols) from min and max of time series
symb_seq_bins = codify(binning, x) #sequence of bin indices -> [1,2,4...]
```
Symbols are represented by integers (index of the bin in which the value falls).
* or for example via ordinal patterns (OP)
```julia
using ComplexityMeasures
op = OrdinalPatterns{3}() #pattern length of 3 (6 possible in total)
symb_seq_op = codify(op, x) #sequence of OP labels -> [1,2,6...]
```

### Get transfer matrix of a symbolic sequence
* by internal method of **StateTransitionNetwork**
```julia
using StateTransitionNetworks
P_bins, Q_bins, ρ_bins =  calculate_transition_matrix(symb_seq_bins)
```
Here generally the returned objects are 
* `P` transfer matrix, where `P[i,j]` is the transition probability from symbol `i -> j` (normalized row-vise: `sum(P[i,:] == 1.0)`)
* `Q` probabilities of transitions, where `Q[i,j]`: transition probability from symbol `i -> j` (normalized sum: `sum(Q) == 1.0)`)
* `ρ`: probabilities of symbol occurrences (from the observed time series)

Transition probabilities from bin to bin will look like:
```
julia> P_bins
4×4 SparseArrays.SparseMatrixCSC{Float64, Int64} with 8 stored entries:
 0.25  0.75   ⋅    ⋅ 
 0.4   0.2   0.2  0.2
  ⋅    1.0    ⋅    ⋅ 
 1.0    ⋅     ⋅    ⋅ 
```
Transition probabilities between ordinal patterns:
```julia
  P_op, Q_op, ρ_op =  calculate_transition_matrix(symb_seq_op)
```
```
julia> P_op
5×5 SparseArrays.SparseMatrixCSC{Float64, Int64} with 7 stored entries:
  ⋅   0.5   ⋅    ⋅   0.5
  ⋅    ⋅   1.0   ⋅    ⋅ 
  ⋅    ⋅    ⋅   1.0   ⋅ 
 0.5  0.5   ⋅    ⋅    ⋅ 
  ⋅    ⋅   1.0   ⋅    ⋅
```

Some functionality from here, for example general transfer matrix construction, was also ported to **ComplexityMeasures**. 
```julia
using ComplexityMeasures
to_op = transferoperator(op,x)
P_op = transfermatrix(to_op)
```
With this, one can access observed ordinal patterns directly
```julia
julia> outcomes(to_op)
  5-element Vector{SVector{3, Int64}}:
   [1, 2, 3]
   [1, 3, 2]
   [3, 2, 1]
   [2, 1, 3]
   [3, 1, 2]
```

### Get dynamic measures using Markov chains: Sinai-Kolmogorov (metric) entropy, Lyapunov measure (trajectory-level varentropy)
```julia
using StateTransitionNetworks
S_KS,_ = sinai_kolmogorov_entropy(P_op)
Λ,_ = lyapunov_measure(P_op)

S_KS, Λ = network_measures(P_op) #together

```
#### Example for a chaotic map
```julia 
using DynamicalSystems
ds = Systems.logistic(r=4.0) 
traj,t = trajectory(ds, 10000; Ttr = 10000)
timeseries = traj[:,1]
op = OrdinalPatterns{5}() #ordinal patterns
to = transferoperator(op,timeseries)
P = transfermatrix(to)
network_measures(P) #KS entropy, Λ for chaotic dynamics (r=4)
 # -> (0.6561118001665418, 0.03863994032705764)

#set to periodic dynamics
ds = Systems.logistic(r=3.5)
traj,t = trajectory(ds, 10000; Ttr = 10000); timeseries = traj[:,1]
to = transferoperator(op,timeseries)
P = transfermatrix(to)
network_measures(P) #KS entropy, Λ for periodic dynamics (r=3.5)
 # -> (0.0, 0.0)

```
### Build transfer matrix incrementally 
```julia
timeseries1 = rand(100) #example time series
timeseries2 = rand(100) #another time series (from a different experiment perhaps)
n_bins = 5
binning = ValueBinning(RectangularBinning(n_bins, true)) #5 bins (symbols) from min and max of time series
symb_seq1 = codify(binning,timeseries1)
symb_seq2 = codify(binning,timeseries2)

#gather all transitions into one object
symb_stats = SymbolStatistics(Int64, n_bins) #expecting symbols of Int64, max number of symbols is n_bins 

#add transitions one-by-one
add_transition!(symb_stats,(1,2)) #add transition 1 -> 2
add_transition!(symb_stats, (2, 1)) #add transition 2 -> 3

#add transitions in one scoop
add_transitions!(symb_stats, symb_seq1) # transitions from timeseries1
add_transitions!(symb_stats,symb_seq2) #transitions from timeseries2

#get transition matrix
P = calculate_transition_matrix(symb_stats)

```

### Using higher-order symbols 
Stringing consecutive symbols together lets us study long-range correlations of symbol transitions. 
```julia
using ComplexityMeasures,StateTransitionNetworks
x = [0.1, 0.2, 0.7, 0.3, 0.8, 0.3, 0.2, 0.8]
binning = ValueBinning(RectangularBinning(3))
symb_seq = codify(binning, x) #first order -> [1,1,3,1,3,1,1,3]
higher_order_symbolics!(symb_seq,2) #turn into pairs of symbols in-place (N-k+1 gets replaced, if k is the order)

julia> symb_seq[1:end-1]
  7-element Vector{Int64}:
    4
    6
   10
    6
   10
    4
    6
```

# References
1. Sándor, B.; Schneider, B.; Lázár, Z.I.; Ercsey-Ravasz M. : [A Novel Measure Inspired by Lyapunov Exponents for the Characterization of Dynamics in State-Transition Networks](https://www.mdpi.com/1099-4300/23/1/103)
2. Sándor, B.; Rusu, A.; Lázár, Z.I.; Ercsey-Ravasz M. : [Measuring dynamical phase transitions in time series](https://doi.org/10.1103/98vk-ws9s)
## Acknowledgements
This work was supported by the grant of the Romanian Ministry of Research, project PN-IV-P2-2.1-TE-2023-1548.
