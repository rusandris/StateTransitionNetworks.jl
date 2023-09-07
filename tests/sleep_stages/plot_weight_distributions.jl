using Plots,LaTeXStrings
using DelimitedFiles
using StateTransitionNetworks
using LinearAlgebra
using StatsBase
using StatsPlots

guidefontsize=20
tickfontsize=15
legendfontsize=15
"figs" in readdir() || mkdir("figs/")


labels = ["NREM2", "SWS", "WAKE", "NREM1"]
colors = [:purple,:green,:red,:orange]



#-----------------netmeasure statistics-----------------------

path_data = "data/netmeasures/"


function plot_EEG_statistics(data;labels,colors,ylabel,tickrotation,alpha)
	pl = violin(1:1:length(data), data, legend=false,color=colors,ylabel= ylabel,
	guidefontsize=guidefontsize,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	dpi=300,
	framestyle=:box,
	ms=1,
	xticks = (1:1:length(data),labels),
	alpha=alpha,
	linewidth=0,
	xrotation=tickrotation)
	
	boxplot!([1:1:length(data);]',data,legend=false,lc=:gray10,color=nothing,lw=2)
	return pl
end

entropies = []
lyapunovs = []

for (i,data_file) in enumerate(readdir(path_data))
	@show path_data
	@show data_file
	rotplane_data = readdlm(path_data*data_file)
	
	push!(entropies,rotplane_data[:,1])
	push!(lyapunovs,rotplane_data[:,2])

end

violin_entr = plot_EEG_statistics(entropies;labels=["NREM2" "SWS" "WAKE" "NREM1"],
	colors=[:purple :green :red :orange],ylabel=L"S_{SK}(\theta)",tickrotation=0,alpha=0.5)

violin_lyap = plot_EEG_statistics(lyapunovs;labels=["NREM2" "SWS" "WAKE" "NREM1"],
	colors=[:purple :green :red :orange],ylabel=L"\Lambda(\theta)",tickrotation=0,alpha=0.5)

plot!(violin_entr,[1,2,3,4],mean.(entropies),mc=[:purple, :green, :red, :orange],st=:scatter,ms=7,ma=0.5,markershape=:utriangle)
plot!(violin_lyap,[1,2,3,4],mean.(lyapunovs),mc=[:purple, :green, :red, :orange],st=:scatter,ms=7,ma=0.5,markershape=:utriangle)

plots_all = plot(violin_entr,violin_lyap,size=(1000,500),layout=(1,2),margin=8Plots.mm)

savefig(plots_all,"figs/sleep_stages_violins.pdf")


#--------------------weight distributions------------------------

path_data = "data/weights_distributions/"

x_axes = [:identity,:identity,:log10,:identity]

distribs = []
distribs_exp = plot()
distribs_power = plot()

labels = ["NREM2", "SWS", "WAKE", "NREM1"]
colors = [:purple,:green,:red,:orange]


for (i,data_file) in enumerate(readdir(path_data))
	@show path_data
	@show data_file
	qweights = readdlm(path_data*data_file)
	h = fit(Histogram,vec(qweights))
	h = normalize(h,mode=:probability)
	bin_centers = h.edges[1][1:end-1] + diff(h.edges[1])/2
	qprob = h.weights
	
	if i == 4
		bin_centers = bin_centers[1:end-4]
		qprob = h.weights[1:end-4]
	end

	
	if i == 2 || i == 4

		plot!(distribs_exp,bin_centers,qprob,
			lw=1,
			ms=6,
			marker=:circle,
			ylabel=L"p(q_{ij})",
			xlabel=L"q_{ij}",
			yaxis=:log10,
			framestyle=:box,
			label=labels[i],
			lc=colors[i],
			mc=colors[i],
			xticks=[0,0.01],
			guidefontsize=guidefontsize,
			tickfontsize=tickfontsize,
			legendfontsize=legendfontsize,reuse=false)
			
	else
		plot!(distribs_power,bin_centers,qprob,
			lw=1,
			ms=6,
			marker=:circle,
			xlabel=L"q_{ij}",
			xaxis=:log10,
			yaxis=:log10,
			lc=colors[i],
			mc=colors[i],
			framestyle=:box,
			label=labels[i],
			guidefontsize=guidefontsize,
			tickfontsize=tickfontsize,
			legendfontsize=legendfontsize,reuse=false)
	end

end


l = @layout [a{0.4h}; b{0.2h}; c{0.2h}; b{0.2h}]
plots_all = plot(distribs_exp,distribs_power,size=(1000,500),layout=(1,2),margin=8Plots.mm)

savefig(plots_all,"figs/sleep_stages_weight_distribs.pdf")
