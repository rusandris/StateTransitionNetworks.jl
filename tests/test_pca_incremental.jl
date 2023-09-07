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


#-----------------params---------------------
data_idxs = 2:4
epoch_length = 386
plane = (1,0.0)
grid = 20
id = 6

data_dir = "./pca_data/pca_selected/"

println("Data reading started from " * data_dir)

#-----------------plotting-------------------

pl_lyaps = plot(
	tickfontsize=12,
	dpi=300,
	c=:red,
	guidefontsize=15,
	legendfontsize=12,
	label=L"\Lambda",
	xlabel=L"epochs ~ added",
	ylabel=L"\Lambda",
	framestyle=:box,
	legend=:outertopright)


pl_entr = plot(
	tickfontsize=12,
	dpi=300,
	c=:blue,
	guidefontsize=15,
	legendfontsize=12,
	label=L"S_{KS}",
	xlabel=L"epoch ~ added",
	ylabel=L"S_{KS}",
	framestyle=:box,
	legend=:outertopright)
	
marker_shapes = [:circle, :rect, :star5, :diamond, :utriangle, :x]
colors = theme_palette(:auto)[1:6]


#-----------------------apply to all files-------------------------------

for (f,data_file) in enumerate(readdir(data_dir))
	@show data_file
	data_filename,ext = splitext(data_file)
	pca_data = readdlm(data_dir*data_file)

	println("Data reading done.")
	println("Shape: ",size(pca_data)," Size in megabytes: ",Base.format_bytes(sizeof(pca_data)))

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



	#discrete timeseries for selected epoch ID
	discrete_timeseries = timeseries_to_common_grid.(epoch_psections[epoch_ids.==id*10], grid, x_min, x_max, y_min, y_max);
	
	for i in 1:length(discrete_timeseries)
		
		stn_added,retcode = add_timeseries(discrete_timeseries[1:i], grid; make_ergodic=true, verbose=false)
		@show nv(stn_added)
		@show retcode
		if retcode == :Success
			P = prob_matrix(stn_added)
			entropy,lyapunov = network_measures(P)
			push!(entropies,entropy)
			push!(lyaps,lyapunov)
			push!(nr_vertices, nv(stn_added))
			push!(densities,ndensity(stn_added))
			
			#plot_stn(stn_added;filename=data_filename*"_"*"stn_added_epoch_id=$(id*10).pdf",nodesize=0.6,nodefillc="orange",linetype="curve",max_edgelinewidth=1)
		end

	end
	
	@show length(lyaps)
	@show length(discrete_timeseries)

	#-----------------plotting measures-------------------

	plot!(pl_lyaps, 1:10:length(lyaps), lyaps[1:10:end], c=:red, label=L"\Lambda"*" : "*data_filename[5:end-4],ms=5,markershape = marker_shapes[f],title="epoch_id=$id")
	plot!(pl_entr, 1:10:length(lyaps), entropies[1:10:end],	c=:blue, label=L"S"*" : "*data_filename[5:end-4],ms=5,markershape = marker_shapes[f],title="epoch_id=$id")
	#plot!(collect(10:10:60), densities,	c=:orange, label=L"\rho")
	#plot!(collect(10:10:60), nr_vertices ./100, c=:green, label=L"nr. vertices/100")
end

savefig(pl_lyaps,"lyaps_incremental_add_epoch_id_$id"*"_test_all_pacients.pdf")
savefig(pl_entr,"entr_incremental_add_epoch_id_$id"*"_test_all_pacients.pdf")
