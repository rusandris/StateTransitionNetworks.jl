module StateTransitionNetworks

import Graphs: DiGraph,add_edge!, inneighbors,outneighbors, nv, ne, weights, is_strongly_connected,strongly_connected_components,edges,degree_histogram,vertices
import StatsBase: sample, mean, var, Weights
using MetaGraphsNext
using GraphPlot,Cairo,Compose
using SparseArrays
using DynamicalSystemsBase: DynamicalSystemsBase,poincaresos,StateSpaceSet,AbstractStateSpaceSet
import LinearAlgebra: eigen, Diagonal, I, nullspace, norm
import DataStructures: OrderedDict
using KrylovKit
export linsolve


include("timeseries_to_grid.jl")
include("transition_matrix.jl")
include("create_stn.jl")
include("iterative_linsolve.jl")
include("network_measures.jl")
include("plot_stn.jl")
include("timeseries_analysis.jl")
include("operations.jl")
include("deprecations.jl")

end
