module StateTransitionNetworks

import Graphs: add_edge!, outneighbors, nv, weights, is_strongly_connected
import SimpleWeightedGraphs: SimpleWeightedDiGraph, get_weight
import StatsBase: sample, mean, var, Weights
import MetaGraphs: MetaDiGraph,set_prop!

include("create_stn.jl")
include("network_measures.jl")

export timeseries_to_grid, create_stn
export walk_statistics, sinai_kolmogorov_entropy,measure_convergence


end
