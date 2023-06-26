using StateTransitionNetworks
using Plots
using LaTeXStrings
using Graphs
using StatsBase
using LinearAlgebra
using DelayEmbeddings
using ChaosTools

include("../plot_functions.jl")

path_data = "./sleep_data/"
#path_data = "./surrogate_data/"
eeg_files = filter(fn -> splitext(fn)[2] == ".dat" || splitext(fn)[2] == ".bin",readdir(path_data))

grid_size = 20
ens = 1000
steps = 10000
plane = (3,0)
var1 = 1
var2 = 4
colors = ["purple","green","red","orange"]
symbol_colors = [:purple,:green,:red,:orange]
labels = Dict("stg2" => "NREM2" ,"stg4" => "SWS","stg5" => "WAKE","stg1" => "NREM1")

plots_all = []


for (i,eeg_file) in enumerate(eeg_files)

	fn,ext = splitext(eeg_file)
	println("\n",fn,"\n")
	
	eeg_data = read_bin(path_data*eeg_file,Float32,6)
	println("Means: ",mean(eeg_data,dims=1))

	
	stn,retcode,entr,lyap = stn_analysis(eeg_data;grid=grid_size,plane=plane,idxs=[var1,var2],return_stn=true)
	@show retcode
		
	@show entr, lyap
	
	density = ndensity(stn)
	@show density
	@show nv(stn)
	
	#plot_degree_distribution!(stn;plot=dd,bins=15,xaxis=:log,yaxis=:log,mc=symbol_colors[i],label=labels[i])

	plot_stn(stn;filename=path_data*fn*".pdf",nodefillc=colors[i],weight_exponent = 1,arrowlengthfrac=0.01,linetype="curve")
	
	
	
	pl_pp = plot_phase_portrait(eeg_data,var1,var2,color=symbol_colors[i])
	push!(plots_all,pl_pp)
	#pl_ts = plot_timeseries(eeg_data,dims=[var1,var2])
	pl_ps = plot_psection(eeg_data;plane = plane,idxs=[var1,var2],color=symbol_colors[i])
	push!(plots_all,pl_ps)

end

	
l = @layout [a{0.4h}; b{0.2h}; c{0.2h}; b{0.2h}]
plots = plot(plots_all...,size=(600,1000),layout=(4,2),left_margin=8Plots.mm,right_margin=8Plots.mm)

savefig(plots,"sleep_stages_plots.pdf")






