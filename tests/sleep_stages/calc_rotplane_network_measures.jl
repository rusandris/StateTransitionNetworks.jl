using Plots
using DelimitedFiles
using StateTransitionNetworks

include("../rotating_plane.jl")

path_data = "sleep_data/"
mkdir("data")
mkdir("data/netmeasures")
mkdir("data/weights_distributions")

#ensemble = 1000
#nsteps = 1000
Δθ = 0.01
grid_size = 20
plane0 = [0.0,0.0,1.0,0.0]

for data_file in readdir(path_data)
	@show data_file
	
	sleep_data = read_bin(path_data*data_file,Float32,6)
	angles,entropies,lyapunovs,PSOS_points,avg_degrees,qweights = 
	rotplane_measures(StateSpaceSet(sleep_data[:,[1,4,3]]);
	grid_size=grid_size,plane0=plane0,θ_min=0,θ_max=π,Δθ=0.01,direction=+1)


	@show qweights[1]
	
	writedlm("data/netmeasures/netmeasures_"*data_file*".txt",hcat(angles,entropies,lyapunovs,PSOS_points,avg_degrees))
	writedlm("data/weights_distributions/qweights_"*data_file*".txt",qweights)
	#savefig(pl,"fig/plot_"*data_file*".png")
end

