using DelimitedFiles
using Plots,LaTeXStrings
#cd("Documents/stn_research/STNResearch/Chaos/revision")
include("../prl_figure_scripts/plotting_params_main.jl")
include("../prl_figure_scripts/plot_functions_chaotic_maps.jl")

#plot logistic
include("plot_EWS_comparison_inset_log.jl")

#plot henon
include("plot_EWS_comparison_inset_henon.jl")

#combine figs 
l = (1,2)
pl = plot(pl_logistic,pl_henon,layout = l,size=(1500,1200)) #size=(1800,1500)
savefig(pl,"EWS_comparison_logistic_henon_inset_denser_od.png")
