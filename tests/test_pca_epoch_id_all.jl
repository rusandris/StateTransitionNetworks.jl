using StateTransitionNetworks
using DelimitedFiles
using StatsBase
using DynamicalSystems
import Graphs: nv,ne
using Plots,LaTeXStrings
using DataFrames,CSV

# cd tests
include("adding_stn_functions.jl")

function translate_to_origin!(epoch; data_idxs=1:3)
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
alldata = DataFrame(Patient=String[],EpochID=Int64[],Entropy=Float64[],Lyapunov=Float64[])

data_dir = "./pca_data/pca_selected/"

println("Data reading started from " * data_dir)

#-----------------plotting measures-IDs-------------------
pl_ids = plot(
	tickfontsize=12,
	dpi=300,
	guidefontsize=15,
	legendfontsize=12,
	xlabel=L"\mathrm{epoch~id}",
	ylabel=L"S,~\Lambda",
	framestyle=:box,
	legend=:outertopright)

	#-----------------plotting Λ-S for epoch IDs-------------------
pl_entrop_lyap = plot(
	tickfontsize=12,
	dpi=300,
	guidefontsize=15,
	legendfontsize=12,
	xlabel=L"\Lambda",
	ylabel=L"S",
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
		
		psection = ChaosTools.poincaresos(Dataset(epoch), plane; direction=+1, idxs=[2,3]);
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
			
			#plot_stn(stn_added;filename=data_filename*"_"*"stn_added_epoch_id=$(id*10).pdf",nodesize=0.6,nodefillc="orange",linetype="curve",max_edgelinewidth=1)
		end
	end

	#-----------------plotting measures-IDs-------------------

	plot!(pl_ids,collect(10:10:60), lyaps, c=:red, label=L"\Lambda"*" : "*data_filename[5:end-4],markershape=marker_shapes[f],ms=5)
	plot!(pl_ids,collect(10:10:60), entropies,	c=:blue, label=L"S"*" : "*data_filename[5:end-4],markershape=marker_shapes[f],ms=5)
	#plot!(collect(10:10:60), densities,	c=:orange, label=L"\rho")
	#plot!(collect(10:10:60), nr_vertices ./100, c=:green, label=L"nr. vertices/100")

	#-----------------plotting Λ-S for epoch IDs-------------------
		
	for i in 1:6
		l = "id=$(i*10)"
		labels = [x[:label] for x in pl_entrop_lyap.series_list]
		l in labels && (l="") 
		scatter!(pl_entrop_lyap,lyaps[i:i],entropies[i:i], label=l,markershape=marker_shapes[f],c=colors[i],ms=5)
	end
	
	#---------------------------create dataframe and add to alldata-----------------------------
	
	patient_data = DataFrame(Patient=fill(data_filename[5:end-4],6),EpochID=10:10:60,Entropy=Float64.(entropies),Lyapunov=Float64.(lyaps))
	@show patient_data
	global alldata = vcat(alldata,patient_data)
	
	

end
savefig(pl_entrop_lyap,"epoch_id_lambda-s_all_pacients.pdf")
savefig(pl_ids,"epoch_id_test_all_pacients.pdf")
CSV.write("pca_stn.csv", alldata)
