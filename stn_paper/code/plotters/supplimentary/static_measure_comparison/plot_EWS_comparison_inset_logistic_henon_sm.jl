using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting_params.jl")
include("../../plot_functions_chaotic_maps.jl")

#plot logistic
include("plot_EWS_comparison_inset_log.jl")

#plot henon
include("plot_EWS_comparison_inset_henon.jl")

fig_dir_name = "figs"
fig_dir = "../../../../" * fig_dir_name * "/" 
mkpath(fig_dir)

#combine figs
l = (1,2)
pl = plot(pl_logistic,pl_henon,layout = l,size=bigfig_size)
savefig(pl,fig_dir*"EWS_comparison_logistic_henon_insets_bigfig.png")
