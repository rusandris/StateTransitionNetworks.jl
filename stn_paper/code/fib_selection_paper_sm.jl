using DelimitedFiles
using CSV,DataFrames
using Plots,LaTeXStrings
using StateTransitionNetworks
using TransitionsInTimeseries
using ComplexityMeasures
using StatsBase
using DataFrames,CSV
cd(@__DIR__)
include("pipeline_functions.jl")
data_dir = "../data/supplimentary/fibrillation_ECG/Long_Term_AF_Database/"
#TODO: mkdir if doesn't exist 
output_dir = "../data/supplimentary/fibrillation_ECG/fibrillation_results/"


#-----------------------------params------------------------
w = 4 #for op
τ = 1 #for op
grid_size = 22 # for binning
#grid_sizes = [20]
grid_edges = [0.2,1.2] # for binning
#grid_edges = Float64[] #let the algo choose the grid minmax
window_size = 400 
lags = [1:10;]

#write out params
writedlm(output_dir*"method_params.txt",[grid_size,w,window_size,grid_edges...])

#4,10,11,12,20 in the paper
samples = ["08","115" ,"39","06"] #filenames near-AF preAF
writedlm(output_dir*"sample_ids.txt",samples)
#-----------------------------read data-------------------------------
#preAF+fibrillation part

rrs_intervals = []
time_spans = []
rrs_datas = []
for i in 1:length(samples)
    rrs_data = readdlm(data_dir*"rr_$(samples[i]).txt")
    rrs = Float64.(rrs_data[2:end,2])
    time_span = [1:length(rrs);] 
    push!(rrs_datas,rrs_data)
    push!(rrs_intervals,rrs)
    push!(time_spans,time_span)
end


#----------------------------------measures and plots--------------------------
#linealphas for different grid grid sizes
#linealphas = range(0.3,1.0;length=length(grid_sizes))
linealphas = [1.0]
plots = []
for i in 1:length(rrs_intervals) 

    #time windows (window stops)
    indicator_window = (width = window_size, stride = 1)
    window_ends = windowmap(last, time_spans[i]; indicator_window...)
    #indicator time series
    measures_grid = windowmap(ts -> early_warning_signals(ts,grid_size,grid_edges,lags;outside_grid=:skip), rrs_intervals[i]; indicator_window...)
    M_grid = stack(measures_grid)'

    measures_OP = windowmap(ts -> early_warning_signals(ts,w,lags;τ=τ), rrs_intervals[i]; indicator_window...)
    M_OP = stack(measures_OP)'

    #write out measures and timeseries data
    #S,Λ,var,ac
    measures_select_grid = M_grid[:,[1,2,7,8]]
    writedlm(output_dir*"measures_grid_$(samples[i])"*"_grid_$grid_size"*"_window_$window_size"*".txt",hcat(window_ends,measures_select_grid))

    measures_select_OP = M_OP[:,[1,2,7,8]]
    writedlm(output_dir*"measures_OP_$(samples[i])"*"_OP_$w"*"_window_$window_size"*".txt",hcat(window_ends,measures_select_OP))
end


