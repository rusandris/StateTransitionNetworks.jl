module StateTransitionNetworks

import Graphs: add_edge!, outneighbors, nv, weights
import SimpleWeightedGraphs: SimpleWeightedDiGraph, get_weight
import StatsBase: sample, mean, var, Weights

include("create_STN.jl")
include("network_measures.jl")

export timeseries_to_grid, create_STN
export random_walk_on_weighted_graph, walk_statistics, sinai_kolmogorov_entropy


end
