module StateTransitionNetworks

import Graphs: SimpleDiGraph,DiGraph,add_edge!, inneighbors,outneighbors, nv, ne, weights, is_strongly_connected,strongly_connected_components,edges,degree_histogram,vertices
import StatsBase: sample, mean, var, Weights
import Base: extrema
using MetaGraphsNext
using GraphPlot,Cairo,Compose
using SparseArrays
using LuxurySparse
using DynamicalSystemsBase: DynamicalSystemsBase,poincaresos,StateSpaceSet,AbstractStateSpaceSet
import DynamicalSystemsBase: minmaxima
import LinearAlgebra: eigen, Diagonal, I, nullspace, norm
import DataStructures: OrderedDict
using KrylovKit
export linsolve


include("timeseries_to_grid.jl")
include("SymbolStatistics.jl")
include("transition_matrix.jl")
include("higher_order_symbolics.jl")

include("iterative_linsolve.jl")
include("network_measures.jl")
include("create_stn.jl")
include("simulate_random_walks.jl")
include("plot_stn.jl")
include("timeseries_analysis.jl")
include("operations.jl")
include("deprecations.jl")

end
