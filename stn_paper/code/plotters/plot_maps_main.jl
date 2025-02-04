using DelimitedFiles
using Plots,LaTeXStrings
using Printf
using DynamicalSystems
#using ComplexityMeasures
include("plotting_params_main.jl")
include("plot_functions_chaotic_maps.jl")

include("plot_logistic_main.jl")
include("plot_henon_main.jl")

l2 = @layout [a{0.5w} b{0.5w}]
pl = plot(pl_logistic,pl_henon,layout=l2,size=colfig_size) #,top_margin=15Plots.mm);

fig_dir_name = "figs"
fig_dir = fig_dir_name * "/" 
!(fig_dir_name in readdir()) && (mkdir(fig_dir))

savefig(pl,fig_dir * "logistic_henon_main.png")
