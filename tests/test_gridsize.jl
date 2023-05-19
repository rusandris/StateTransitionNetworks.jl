using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using LaTeXStrings
using StatsBase


function ndensity(stn)
	nr_vertices = nv(stn)
	ne(stn)/(nr_vertices*(nr_vertices-1))
end


T = 5000
rho_vals = [180.1, 180.7, 180.78]
plane = (1,15.0)
ds = Systems.lorenz();
grid_sizes = 5:150
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
			push!(entropies,entropy)
			push!(lyaps,lyapunov)
			push!(degrees, mean(outdegree(stn)))
			push!(densities,ndensity(stn))
		end
		
	end
	
	
	pl_rho = plot(lyaps,
	tickfontsize=12,
	dpi=300,
	c=:red,
	guidefontsize=15,
	legendfontsize=12,
	label=L"\Lambda",
	xlabel=L"grid ~ size",
	ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
	framestyle=:box)


	plot!(entropies,
		tickfontsize=12,
		dpi=300,
		c=:blue,
		guidefontsize=15,
		legendfontsize=12,
		label=L"S_{KS}",
		xlabel=L"grid ~ size",
		ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
		framestyle=:box)
		
	plot!(densities,
		tickfontsize=12,
		dpi=300,
		c=:orange,
		guidefontsize=15,
		legendfontsize=12,
		label=L"\rho : density",
		xlabel=L"grid ~ size",
		ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
		framestyle=:box)


	
	plot!(degrees ./10,
		tickfontsize=12,
		dpi=300,
		c=:green,
		guidefontsize=15,
		legendfontsize=12,
		label=L"<k_{out}>/10",
		xlabel=L"grid ~ size",
		ylabel=L"nv,~\rho,~S_{KS},~\Lambda",
		framestyle=:box,
		legend=:outertopright)
		

		
	    savefig(pl_rho, "grid_size_effect_T_$T"*"_rho_$(rho).svg")
	
end
