using DelimitedFiles
using CSV,DataFrames
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting_params.jl")

include("plot_noise_results_logistic.jl")
include("plot_noise_results_henon.jl")

fig_dir_name = "figs"
fig_dir = "../../../../" * fig_dir_name * "/" 
!(fig_dir_name in readdir("../../../../")) && (mkdir(fig_dir))

savefig(pl_log,fig_dir*"logistic_noise_results.png")
savefig(pl_henon,fig_dir*"/henon_noise_results.png")
