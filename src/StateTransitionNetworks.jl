module StateTransitionNetworks

import Graphs: add_edge!, outneighbors, nv, weights, is_strongly_connected,edges
import SimpleWeightedGraphs: SimpleWeightedDiGraph, get_weight
import StatsBase: sample, mean, var, Weights
import MetaGraphs: MetaDiGraph,set_prop!,get_prop,set_indexing_prop!
using GraphPlot,Cairo,Compose

include("create_stn.jl")
include("network_measures.jl")
include("plot_stn.jl")

export timeseries_to_grid, create_stn
export walk_statistics, sinai_kolmogorov_entropy,measure_convergence
export plot_stn


end
