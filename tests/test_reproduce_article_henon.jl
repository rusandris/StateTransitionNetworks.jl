using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using LaTeXStrings
using StatsBase
using DelimitedFiles


T = 30000
Ttr = 1000
a_vals = [1:0.01:1.4;]
ds = Systems.henon();
grid_size = 20
Δt = 0.01
traj_ensemble = 1
rw_ensemble = 100 
nr_steps = 100


analytic_lyapunovs = []
analytic_entropies = []
num_lyapunovs = []
num_entropies = []


for a in a_vals
	@show a

	set_parameter!(ds,1,a)
	analytic_entropy = 0.0
	analytic_lyapunov = 0.0
	num_entropy = 0.0
	num_lyapunov = 0.0
	
	
	for i in 1:traj_ensemble
		@show i
		x0 = rand() .- 0.5 
		y0 = rand() * 0.4 .- 0.2 
		timeseries = trajectory(ds, T,[x0,y0]; Ttr=Ttr);
		d_traj, v_names = timeseries_to_grid(timeseries, grid_size);
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

writedlm("data_reproduced_netmeasures_henon_a_T_$T"*"_Ttr_$Ttr"*"_traj_ens_$traj_ensemble"*"_rwens_$rw_ensemble"*"_nsteps_$nr_steps"*"_grid_$grid_size"*"standardpsos.txt",hcat(a_vals,analytic_entropies,analytic_lyapunovs,num_entropies,num_lyapunovs))

