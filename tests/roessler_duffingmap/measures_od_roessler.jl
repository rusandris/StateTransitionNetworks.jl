using StateTransitionNetworks
using DynamicalSystems
using Graphs
using StatsBase
using DelimitedFiles

T = 5000
nr_points = 5000
Ttr = 500
b_vals = [0.2:0.01:1.0;]

#stn_colors = Dict(zip([180.10,180.70,180.78,181.10],[colorant"orange",colorant"red",colorant"green",colorant"blue"]))



plane = (2,0.0)
ds = PredefinedDynamicalSystems.roessler();
pmap = PoincareMap(ds,plane;rootkw=(xrtol=1e-10, atol=1e-10), direction=1)
grid_size = 20
Δt = 0.01
rw_ensemble = 100 
nr_steps = 10000


#---------------orbit diagram-----------------------

println("Orbit diagram calculation starting...")
od = orbitdiagram(pmap, 1, 2, b_vals; n=nr_points, Ttr=Ttr, show_progress=true)
writedlm("data/orbit_diagram_roessler_b_saved_z_T_$T"*"_Ttr_$Ttr.txt",od)
println("Done.")

#---------------normal lyapunovs--------------------

println("Lyapunov exponent calculation starting...")
lyapunov_exponents = zeros(length(b_vals))
for (i,b) in enumerate(b_vals)
	@show b

	set_parameter!(ds,2,b)
	
	lyap_exp = lyapunov(ds,T;Ttr = Ttr)
	lyapunov_exponents[i] = lyap_exp

end

writedlm("data/lyapunov_exponent_roessler_b_T_$T"*"_Ttr_$Ttr.txt",hcat(b_vals,lyapunov_exponents))
println("Done.")



#-----------------calculate network measures----------------------


analytic_lyapunovs = []
analytic_entropies = []
num_lyapunovs = []
num_entropies = []


for b in b_vals
	@show b

	set_parameter!(ds,2,b)

	timeseries,  = trajectory(ds, T; Δt=Δt, Ttr=Ttr);
	stn, retcode = create_stn(timeseries, grid_size, plane, [1,3];make_ergodic=true,verbose=true,direction=+1);
	
	if retcode == :Success
		
		P = prob_matrix(stn)
		
	
		analytic_measures = network_measures(P)
		num_measures = network_measures(stn,rw_ensemble,nr_steps)
		
		push!(analytic_entropies,analytic_measures[1])
		push!(analytic_lyapunovs, analytic_measures[2])
		push!(num_entropies,num_measures[1])
		push!(num_lyapunovs, num_measures[2])
	else
		@warn "STN creation wasn't successful!"
		
		push!(analytic_entropies,NaN)
		push!(analytic_lyapunovs, NaN)
		push!(num_entropies,NaN)
		push!(num_lyapunovs, NaN)
	end
	

		
end

writedlm("data/data_netmeasures_roessler_b_T_$T"*"_Ttr_$Ttr"*"_traj_ens"*"_rwens_$rw_ensemble"*"_nsteps_$nr_steps"*"_grid_$grid_size"*"standardpsos.txt",hcat(b_vals,analytic_entropies,analytic_lyapunovs,num_entropies,num_lyapunovs))



