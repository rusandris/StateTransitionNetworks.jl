using StateTransitionNetworks
using DynamicalSystems
using StatsBase
using DelimitedFiles
using Printf
cd(@__DIR__)

include("functions_chaotic_maps.jl")

data_dir = "../data/supplimentary/roessler_data/" 
mkpath(data_dir)

T = 10^4 #length of time series (nr of poinc points) for measures
nr_points = 100 #nr of poincare points on od
Ttr = 10^4 #10000
b_vals = [0.35:0.0001:0.4;]
b_start = b_vals[begin]
b_stop = b_vals[end]

T_string::String = @sprintf "%.E" T
Ttr_string::String = @sprintf "%.E" Ttr

plane = (2,0.0)
plane_string::String = "$(plane[1])"*"_"*"$(plane[2])"

ds = Systems.roessler();
pmap = PoincareMap(ds,plane;rootkw=(xrtol=1e-10, atol=1e-10), direction=1) #becomes discrete system
pmap = ProjectedDynamicalSystem(pmap,[1,3],[plane[2]]) #take out the the reduced dimension (3d->2d)
grid_size = 2^5
grid_edges::Vector{Float64} = [-10.0,0.0,-2.0,0.5]
Δt = 0.01
rw_ensemble = 1000 
nr_steps = 10000

orders::Vector{Int64} = [1,2,3,4] #[1,4,8,12]
const n::Int64 = 100 
const N_steps::Int64 = 500
const ensemble::Int64 = 100 #1000
const ϵ::Float64 = 1e-3 #1e-5

out_file_entropy = data_dir*"roessler_entropies"  * "_n_poinc_$nr_points" * "plane_$plane_string" * "_b_0.3" * "_grid_$grid_size" * "_param_$b_start" * "_$b_stop" * ".txt"
out_file_lambda = data_dir*"roessler_lambdas" *  "_n_poinc_$nr_points" * "plane_$plane_string" * "_b_0.3" * "_grid_$grid_size" * "param_$b_start" * "_$b_stop" * ".txt"


#---------------orbit diagram-----------------------

println("Orbit diagram calculation starting...")
od = orbitdiagram(pmap, 1, 2, b_vals; n=nr_points, Ttr=Ttr, show_progress=true)
writedlm(data_dir*"od_roessler_b_saved_z_T_$T"*"_Ttr_$Ttr"*"_b_$(b_start)_$b_stop"*".txt",stack(od)')
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

writedlm(data_dir*"lyapexps_roessler_b_T_$T"*"_Ttr_$Ttr"*"_b_$(b_start)_$b_stop"*".txt",hcat(b_vals,lyapunov_exponents))
println("Done.")



#-----------------calculate measures----------------------


analytic_lyapunovs = []
analytic_entropies = []
num_lyapunovs = []
num_entropies = []

#--------------------------calc_measures------------------------

writedlm(data_dir * "orders" * ".txt",orders)
writedlm(data_dir * "roessler_p_values" * "_nr_param_$(length(b_vals))"*".txt",b_vals)

@info "Calculating measures...."
flush(stderr)

calc_measures(pmap,b_vals,orders;
	sts_eltype=Int128,
	out_file_entropy=out_file_entropy,
	out_file_lambda=out_file_lambda,
	param_index=2,
	grid_size=grid_size,
	grid_edges=grid_edges,
	u0=initial_state(ds),
	T=T,
	Ttr=Ttr,
	ensemble = ensemble,
	N_steps = N_steps,
	ϵ=ϵ)
@info "Measure calculation done and saved."
flush(stderr)




writedlm("data/data_netmeasures_roessler_b_T_$T"*"_Ttr_$Ttr"*"_traj_ens"*"_rwens_$rw_ensemble"*"_nsteps_$nr_steps"*"_grid_$grid_size"*"standardpsos" *"_b_$(b_start)_$b_stop"*".txt",hcat(b_vals,analytic_entropies,analytic_lyapunovs,num_entropies,num_lyapunovs))



