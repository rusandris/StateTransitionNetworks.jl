using StateTransitionNetworks
using DynamicalSystems
using Graphs
using StatsBase
using DelimitedFiles


T = 30000
Ttr = 1000
a_vals = [2.505:0.0005:2.51;]

#stn_colors = Dict(zip([180.10,180.70,180.78,181.10],[colorant"orange",colorant"red",colorant"green",colorant"blue"]))


duffingmap(u,p,n) = SVector(u[2],-p[2]*u[1] +p[1]*u[2] - u[2]^3) 
p = [2.3,0.1]
u0 = [0.5,0.2]

ds = DeterministicIteratedMap(duffingmap,u0,p);
grid_size = 20
rw_ensemble = 100 
nr_steps = 10000



#---------------orbit diagram-----------------------

println("Orbit diagram calculation starting...")
od = orbitdiagram(ds, 1, 1, a_vals; n=1000, Ttr=Ttr)
writedlm("data/orbit_diagram_duffingmap_a_saved_x_T_$T"*"_Ttr_$Ttr.txt",od)
println("Done.")

#---------------normal lyapunovs--------------------

println("Lyapunov exponent calculation starting...")
lyapunov_exponents = zeros(length(a_vals))
for (i,a) in enumerate(a_vals)
	@show a

	set_parameter!(ds,1,a)
	
	lyap_exp = lyapunov(ds,T;Ttr = Ttr)
	lyapunov_exponents[i] = lyap_exp

end

writedlm("data/lyapunov_exponent_duffingmap_a_T_$T"*"_Ttr_$Ttr.txt",hcat(a_vals,lyapunov_exponents))
println("Done.")

#----------------calculate network measures-----------------


analytic_lyapunovs = []
analytic_entropies = []
num_lyapunovs = []
num_entropies = []


for a in a_vals
	@show a

	set_parameter!(ds,1,a)
	
	timeseries,  = trajectory(ds, T;Ttr=Ttr);
	stn, retcode = create_stn(timeseries,grid_size;make_ergodic=true,verbose=true);
	
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

writedlm("data/data_netmeasures_duffingmap_b_T_$T"*"_Ttr_$Ttr"*"_rwens_$rw_ensemble"*"_nsteps_$nr_steps"*"_grid_$grid_size"*".txt",hcat(a_vals,analytic_entropies,analytic_lyapunovs,num_entropies,num_lyapunovs))





