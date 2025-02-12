using DelimitedFiles
using Plots,LaTeXStrings
using Printf
using DynamicalSystems
#using ComplexityMeasures
cd(@__DIR__)
include("../../plotting_params.jl")
include("../../plot_functions_chaotic_maps.jl")

#include separate plotting scripts, merge everything here

include("plot_logistic_sm.jl")
include("plot_henon_sm.jl")

l2 = @layout [a{0.5w} b{0.5w}]
pl = plot(pl_logistic,pl_henon,layout=l2,size=bigfig_size);

fig_dir_name = "figs"
fig_dir = "../../../../" * fig_dir_name * "/" 
!(fig_dir_name in readdir("../../../../")) && (mkdir(fig_dir))
savefig(pl,fig_dir * "logistic_henon_sm.png")

