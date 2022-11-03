using DynamicalSystems
using StateTransitionNetworks
using Plots
using LaTeXStrings
using OrdinaryDiffEq	

ds = Systems.lorenz()
tol = 1e-8
params = [180.7,180.1,180.78,181.1]
ensemble = 10000
N_max = 5000
plane = (1,15.0)
T = 5000
Ttr = 500
steps = 1:N_max
colors = [:red,:orange,:green,:blue]
grid_size = 20
#diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true)


pl = plot()

for (i,ρ) in enumerate(params)
	set_parameter!(ds,2,ρ)
	@show ρ
	psection = ChaosTools.poincaresos(ds, plane, T; Ttr=Ttr, direction=+1,rootkw = (xrtol = tol, atol = tol));
	traj = psection[:,2:end] 
	traj_grid, vertex_names = timeseries_to_grid(traj,grid_size) 
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
		label="$ρ",
		legendfontsize=15,
		legend=:best,
		dpi=300)
	end
end

savefig(pl,"lorenz_convergence_ens$ensemble"*"_N$N_max"*"_T$T"*"_Ttr$Ttr")
