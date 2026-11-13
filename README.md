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
  symb_seq_bins = codify(binning, x) #symbolic time series with binning
```
Symbols are represented by integers (index of the bin in which the value falls).
* or for example via ordinal patterns (OP)
```julia
  using ComplexityMeasures
  op = OrdinalPatterns{3}() #pattern length of 3 (6 possible in total)
  symb_seq_op = codify(op, x) #symbolic time series with OP
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
  to_op = transferoperator(op,symb_seq_op)
  P_op = transfermatrix(to_op)
```
### Get dynamic measures using Markov chains: Sinai-Kolmogorov (metric) entropy, Lyapunov measure (trajectory-level varentropy)
```julia
  using StateTransitionNetworks
  S_KS,_ = sinai_kolmogorov_entropy(P_op)
  Λ,_ = lyapunov_measure(P_op)

  S_KS, Λ = network_measures(P_op) #together

```


**TODO** : 
  - add logistic/henon example for KS,Lambda
  - add incremental construction example
  - add reference to PRR publication
  - show symbols?

# References
Sándor, B.; Schneider, B.; Lázár, Z.I.; Ercsey-Ravasz M. : [A Novel Measure Inspired by Lyapunov Exponents for the Characterization of Dynamics in State-Transition Networks](https://www.mdpi.com/1099-4300/23/1/103)
