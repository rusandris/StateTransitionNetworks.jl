using StateTransitionNetworks
using DelimitedFiles
using StatsBase
using DynamicalSystems
import Graphs: nv,ne
using Plots,LaTeXStrings

# cd tests
include("adding_stn_functions.jl")

function translate_to_origin!(epoch; data_save_idxs=1:3)
	for i in data_idxs
		epoch[:,i] = epoch[:,i] .- mean(epoch[:,i])
	end
end


function ndensity(stn)
	nr_vertices = nv(stn)
	ne(stn)/(nr_vertices*(nr_vertices-1))
end


#read data
pca_data = readdlm("pca_data/sub-AH316_pca.txt")

#-----------------params---------------------
data_idxs = 2:4
epoch_length = 386
plane = (1,0.0)
grid = 20


#------------------get psos for every epoch------------------------------

epoch_psections = []


for epoch_start_index in 1:epoch_length:(size(pca_data)[1]-epoch_length)
	epoch = pca_data[epoch_start_index:epoch_start_index+epoch_length-1,data_idxs]
	#@show epoch_start_index
	#@show size(epoch)[1]
	
	translate_to_origin!(epoch)
	
	psection = DynamicalSystemsBase.poincaresos(Dataset(epoch), plane; direction=+1, save_idxs=[2,3]);
	push!(epoch_psections,psection)
	
end

#common grid extrema
x_min, x_max, y_min, y_max = get_grid_edges(epoch_psections)

#list of discrete timeseries for every epoch
discrete_timeseries = timeseries_to_common_grid.(epoch_psections, grid, x_min, x_max, y_min, y_max);



lyaps = []
entropies = []
densities = []
nr_vertices = []


#----------------------------create stns by adding epochs---------------------------------------

for i in 1:length(discrete_timeseries)
	stn_added,retcode = add_timeseries(discrete_timeseries[1:i],grid;make_ergodic=true,verbose=false)
	@show nv(stn_added)
	@show retcode
	P = prob_matrix(stn_added)
	if retcode == :Success
		entropy,lyapunov = network_measures(P)
		push!(entropies,entropy)
		push!(lyaps,lyapunov)
		push!(nr_vertices, nv(stn_added))
		push!(densities,ndensity(stn_added))
		
		if i == length(discrete_timeseries)
			plot_stn(stn_added;filename="stn_added.pdf",nodesize=0.6,nodefillc="orange",linetype="curve",max_edgelinewidth=1)
		end
		
	end
end


#-----------------plotting-------------------

pl = plot(lyaps,
	tickfontsize=12,
	dpi=300,
	c=:red,
	guidefontsize=15,
	legendfontsize=12,
	label=L"\Lambda",
	xlabel=L"epochs ~ added",
	ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
	framestyle=:box)


plot!(entropies,
	tickfontsize=12,
	dpi=300,
	c=:blue,
	guidefontsize=15,
	legendfontsize=12,
	label=L"S_{KS}",
	xlabel=L"epoch ~ added",
	ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
	framestyle=:box)
	
plot!(densities,
	tickfontsize=12,
	dpi=300,
	c=:orange,
	guidefontsize=15,
	legendfontsize=12,
	label=L"\rho : density",
	xlabel=L"epochs ~ added",
	ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
	framestyle=:box)


plot!(nr_vertices ./100,
	tickfontsize=12,
	dpi=300,
	c=:green,
	guidefontsize=15,
	legendfontsize=12,
	label=L"nr. vertices/100",
	xlabel=L"epochs ~ added",
	ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
	framestyle=:box)

savefig(pl,"add_stn_epoch_test.pdf")



#-----------------plot epoch-stns-------------------


for i in 1:10
	stn_added,retcode = add_timeseries(discrete_timeseries[i:i],grid;make_ergodic=true,verbose=true)
	@show nv(stn_added)
	@show retcode
	
	if retcode == :Success
		P = prob_matrix(stn_added)
		entropy,lyapunov = network_measures(P)
		plot_stn(stn_added;filename="stn_pca_epoch$(i).pdf",nodesize=0.6,nodefillc="orange",linetype="curve",max_edgelinewidth=1)
	end
end

