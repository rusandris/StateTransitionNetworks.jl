using DelimitedFiles
using Plots,LaTeXStrings
using Printf
using DynamicalSystems
#using ComplexityMeasures
cd(@__DIR__)
include("plotting_params.jl")
include("plot_functions_chaotic_maps.jl")

#main fig
include("main/plot_maps_main.jl")

#supplimentary
include("supplimentary/fibrillation_ECG/plot_fib_selection_paper_sm.jl")
include("supplimentary/noisey_timeseries/plot_noise_results.jl")
include("supplimentary/sliding_param_measures/plot_slideparam.jl")
include("supplimentary/static_measure_comparison/plot_EWS_comparison_inset_logistic_henon.jl")
include("supplimentary/wider_range_measures/plot_maps_sm.jl")
