using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using LaTeXStrings
using StatsBase
using DelimitedFiles


T = 5000
Ttr = 1000
rho_vals = [180:0.01:182.0;]
plane = (1,15.0)
ds = PredefinedDynamicalSystems.lorenz();
grid_size = 20
Δt = 0.01
traj_ensemble = 1
rw_ensemble = 100 
nr_steps = 100


analytic_lyapunovs = []
analytic_entropies = []
num_lyapunovs = []
num_entropies = []


for rho in rho_vals
	@show rho

	set_parameter!(ds,2,rho)
	analytic_entropy = 0.0
	analytic_lyapunov = 0.0
	num_entropy = 0.0
	num_lyapunov = 0.0
	
	
	for i in 1:traj_ensemble
		@show i
		timeseries,  = trajectory(ds, T, rand(3)*20 .- 10;Δt=Δt, Ttr=Ttr);
		psection = DynamicalSystemsBase.poincaresos(timeseries, plane; direction=+1, save_idxs=[2,3]);
		d_traj, v_names = timeseries_to_grid(psection, grid_size);
		stn, retcode = create_stn(d_traj, v_names;make_ergodic=true,verbose=true);
		
		P = prob_matrix(stn)
		if retcode == :Success
			analytic_measures = network_measures(P)
			num_measures = network_measures(stn,rw_ensemble,nr_steps)
			
			analytic_entropy += analytic_measures[1]
			analytic_lyapunov += analytic_measures[2]
			num_entropy += num_measures[1]
			num_lyapunov += num_measures[2]
		else
			@warn "STN creation wasn't successful!"
		end
		
	end
		
	push!(analytic_entropies,analytic_entropy/traj_ensemble)
	push!(analytic_lyapunovs,analytic_lyapunov/traj_ensemble)
	push!(num_entropies,num_entropy/traj_ensemble)
	push!(num_lyapunovs,num_lyapunov/traj_ensemble)
		
end

writedlm("data_reproduced_netmeasures_lorenz_rho_T_$T"*"_Ttr_$Ttr"*"_traj_ens_$traj_ensemble"*"_rwens_$rw_ensemble"*"_nsteps_$nr_steps"*"_grid_$grid_size"*"standardpsos.txt",hcat(rho_vals,analytic_entropies,analytic_lyapunovs,num_entropies,num_lyapunovs))

