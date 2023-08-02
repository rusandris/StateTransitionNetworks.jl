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
epoch_ids = []

for epoch_start_index in 1:epoch_length:(size(pca_data)[1]-epoch_length)
	epoch = pca_data[epoch_start_index:epoch_start_index+epoch_length-1,data_idxs]
	epoch_id = pca_data[epoch_start_index:epoch_start_index+epoch_length-1,1]
	#@show epoch_start_index
	#@show size(epoch)[1]
	
	translate_to_origin!(epoch)
	
	psection = DynamicalSystemsBase.poincaresos(Dataset(epoch), plane; direction=+1, save_idxs=[2,3]);
	push!(epoch_psections,psection)
	push!(epoch_ids, Int(unique(epoch_id)[1]))
	
end

#common grid extrema
x_min, x_max, y_min, y_max = get_grid_edges(epoch_psections)

# collect data based on their epoch IDs
lyaps = []
entropies = []
densities = []
nr_vertices = []
for id in 1:6
	#list of discrete timeseries for every epoch ID
	discrete_timeseries = timeseries_to_common_grid.(epoch_psections[epoch_ids.==id*10], grid, x_min, x_max, y_min, y_max);
	stn_added,retcode = add_timeseries(discrete_timeseries, grid; make_ergodic=true, verbose=false)
	@show nv(stn_added)
	@show retcode
	if retcode == :Success
		P = prob_matrix(stn_added)
		entropy,lyapunov = network_measures(P)
		push!(entropies,entropy)
		push!(lyaps,lyapunov)
		push!(nr_vertices, nv(stn_added))
		push!(densities,ndensity(stn_added))
		
		plot_stn(stn_added;filename="stn_added_epoch_id=$(id*10).pdf",nodesize=0.6,nodefillc="orange",linetype="curve",max_edgelinewidth=1)
	end
end

#-----------------plotting measures-IDs-------------------
pl = plot(
	tickfontsize=12,
	dpi=300,
	guidefontsize=15,
	legendfontsize=12,
	xlabel=L"\mathrm{epoch~id}",
	ylabel=L"S,~\Lambda",
	framestyle=:box)
plot!(collect(10:10:60), lyaps, c=:red, label=L"\Lambda")
plot!(collect(10:10:60), entropies,	c=:blue, label=L"S")
#plot!(collect(10:10:60), densities,	c=:orange, label=L"\rho")
#plot!(collect(10:10:60), nr_vertices ./100, c=:green, label=L"nr. vertices/100")
savefig(pl,"epoch_id_test.pdf")

#-----------------plotting Λ-S for epoch IDs-------------------
pl = plot(
	tickfontsize=12,
	dpi=300,
	guidefontsize=15,
	legendfontsize=12,
	xlabel=L"\Lambda",
	ylabel=L"S",
	framestyle=:box)
for i in 1:6
	scatter!(pl, lyaps[i:i], entropies[i:i], label="id=$(i*10)")
end

savefig(pl,"epoch_id_lambda-s.pdf")
