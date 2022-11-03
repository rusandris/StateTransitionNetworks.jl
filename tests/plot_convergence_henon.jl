using DynamicalSystems
using StateTransitionNetworks
using Plots
using LaTeXStrings

ds = Systems.henon()
params = [1.2,1.2265,1.24,1.27]
ensemble = 10000
N_max = 5000
T = 10000
Ttr = 1000
steps = 1:N_max
colors = [:orange,:red,:blue,:green]

pl = plot()

for (i,a) in enumerate(params)
	set_parameter!(ds,1,a)
	@show a
	traj = trajectory(ds,T;Ttr = Ttr)
	traj_grid, vertex_names = timeseries_to_grid(traj,20)
	stn,retcode = create_stn(traj_grid,vertex_names)
	if retcode ==:Success
		entr,lyaps = measure_convergence(stn,ensemble,N_max)
				
		plot!(steps,lyaps,
		lw=2,
		lc = colors[i],
		xlabel="t",
		ylabel=L"$\Delta L/t$",
		xguidefontsize=18,
		yguidefontsize=18,
		tickfontsize=10,
		label="$a",
		legendfontsize=15,
		legend=:best,
		dpi=300)
	end
end

savefig(pl,"henon_convergence_ens$ensemble"*"_N$N_max"*"_T$T"*"_Ttr$Ttr")
