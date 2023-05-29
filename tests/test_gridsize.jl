using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using LaTeXStrings
using StatsBase
using DelimitedFiles


function ndensity(stn)
	nr_vertices = nv(stn)
	ne(stn)/(nr_vertices*(nr_vertices-1))
end


T = 50
rho_vals = [180.1, 180.7, 180.78]
plane = (1,15.0)
ds = Systems.lorenz();
grid_sizes = 5:20
gridmax = grid_sizes[end]
Δt = 0.01



for rho in rho_vals
	@show rho
	lyaps = []
	entropies = []
	degrees = []
	densities = []

	set_parameter!(ds,2,rho)
	timeseries = trajectory(ds, T; Δt=Δt, Ttr=500);
	psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
	@show length(psection)
	
	
	for grid_size in grid_sizes
		@show grid_size
		d_traj, v_names = timeseries_to_grid(psection, grid_size);
		stn, retcode = create_stn(d_traj, v_names;make_ergodic=true,verbose=true);
		
		P = prob_matrix(stn)
		if retcode == :Success
			entropy,lyapunov = network_measures(P)
			#entropy,lyapunov = network_measures(stn,1000,1000)
			push!(entropies,entropy)
			push!(lyaps,lyapunov)
			push!(degrees, mean(outdegree(stn)))
			push!(densities,ndensity(stn))
		end
		
	end
	
	
	pl_rho = plot(lyaps,
	tickfontsize=15,
	dpi=300,
	c=:red,
	guidefontsize=18,
	legendfontsize=15,
	label=L"\Lambda",
	xlabel=L"grid ~ size",
	ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
	framestyle=:box)


	plot!(entropies,
		tickfontsize=15,
		dpi=300,
		c=:blue,
		guidefontsize=18,
		legendfontsize=15,
		label=L"S_{KS}",
		xlabel=L"grid ~ size",
		ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
		framestyle=:box)
		
	plot!(densities,
		tickfontsize=15,
		dpi=300,
		c=:orange,
		guidefontsize=18,
		legendfontsize=15,
		label=L"\rho : density",
		xlabel=L"grid ~ size",
		ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
		framestyle=:box)


	
	plot!(degrees ./10,
		tickfontsize=15,
		dpi=300,
		c=:green,
		guidefontsize=18,
		legendfontsize=15,
		label=L"<k_{out}>/10",
		xlabel=L"grid ~ size",
		ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
		framestyle=:box,
		legend=:outertopright)
		
	plot!(title=L"\rho = "*"$rho",titlefontsize=15)
	
	writedlm("data_grid_size_effect_T_$T"*"_rho_$(rho)"*"_gridmax=$gridmax"*".txt",hcat(lyaps,entropies,densities,degrees))	
	savefig(pl_rho, "grid_size_effect_T_$T"*"_rho_$(rho)"*"_gridmax=$gridmax"*".svg")
	
end
