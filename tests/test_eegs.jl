using StateTransitionNetworks
using Plots
using LaTeXStrings
using Graphs
using StatsBase
using LinearAlgebra
using DelayEmbeddings
using ChaosTools

include("plot_functions.jl")

eeg_files = filter(fn -> splitext(fn)[2] == ".dat" || splitext(fn)[2] == ".bin",readdir())

grid = 20
ens = 1000
steps = 10000
plane = (1,0)
var1 = 2
var2 = 4
colors = ["purple","green","red","orange"]
symbol_colors = [:purple,:green,:red,:orange]
labels = ["sgt2","sgt4","sgt5","sgt1"]

dd = plot()
for (i,eeg_file) in enumerate(eeg_files)
	eeg_data = read_bin(eeg_file,Float32,6)
	fn,ext = splitext(eeg_file)
	println(fn)
	
	stn,retcode = stn_analysis(eeg_data;grid=grid,plane=plane,idxs=[var1,var2],return_stn=true)
	@show retcode

	
	entr, lyap = network_measures(stn,ens,steps)
	@show entr, lyap
	
	density = ndensity(stn)
	@show density
	
	
	plot_degree_distribution!(stn;plot=dd,bins=15,xaxis=:log,yaxis=:log,mc=symbol_colors[i],label=labels[i])
	
	plot_stn(stn;filename=fn*".pdf",nodefillc=colors[i])
	
	plot_phase_portrait(eeg_data,var1,var2;filename = "pp_vars$var1"*"$var2" * fn * ".png")
	plot_timeseries(eeg_data,filename="ts_vars$var1"*"$var2" * fn * ".png",dims=[var1,var2])
	plot_psection(eeg_data;plane = plane,idxs=[var1,var2],filename="psec_vars$var1"*"$var2" * fn * ".png")
	
end

savefig(dd,"eeg_degree_distributions.png")





