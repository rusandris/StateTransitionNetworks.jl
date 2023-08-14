using StateTransitionNetworks
using DynamicalSystems
using Graphs
using StatsBase
using DelimitedFiles
using Plots,LaTeXStrings

#create data and fig dirs
"data" in readdir() || mkdir("data/")
"figs" in readdir() || mkdir("figs/")

#run scripts

include("measures_od_duffingmap.jl")
include("measures_od_roessler.jl")
include("plot_measures_duffingmap.jl")
include("plot_measures_roessler.jl")
include("plot_stns_trajs_roessler.jl")
include("plot_stns_trajs_duffingmap.jl")
include("plot_walklength_distribution.jl")
include("test_convergence_roessler.jl")
include("test_rotating_plane_roessler.jl")
include("test_rotating_plane_roessler.jl")
include("roessler_psection3D.jl")
include("rotplane_measures_lorenz.jl")
include("plot_rotplane_measures_lorenz.jl")




