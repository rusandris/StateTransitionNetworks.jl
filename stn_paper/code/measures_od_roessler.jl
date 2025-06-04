using StateTransitionNetworks
using DynamicalSystems
using StatsBase
using OrdinaryDiffEq
using DelimitedFiles
using Printf
cd(@__DIR__)

include("functions_chaotic_maps.jl")

data_dir = "../data/supplimentary/roessler_data/" 
mkpath(data_dir)

T = 5*10^5 #length of time series (nr of poinc points) for measures
nr_points = 100 #nr of poincare points on od
Ttr = 10^4 #10000
b_vals = [0.35:0.0001:0.4;]
#b_vals = [0.368:0.00005:0.4;]
b_start = b_vals[begin]
b_stop = b_vals[end]

T_string::String = @sprintf "%.E" T
Ttr_string::String = @sprintf "%.E" Ttr

plane = (2,0.0)
plane_string::String = "$(plane[1])"*"_"*"$(plane[2])"
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-1,maxiters=1e15)

#roessler system
@inbounds function roessler_rule(u, p, t)
    du1 = -u[2]-u[3]
    du2 = u[1] + p[1]*u[2]
    du3 = p[2] + u[3]*(u[1] - p[3])
    return SVector{3}(du1, du2, du3)
end

u0=[1.0, -2.0, 0.1]
p = [0.2,0.2,5.7] #a = 0.2, b = 0.2, c = 5.7
ds = CoupledODEs(roessler_rule, u0, p;diffeq=diffeq)


pmap = PoincareMap(ds,plane; rootkw=(xrtol=1e-10, atol=1e-10), direction=1) #becomes discrete system
pmap = ProjectedDynamicalSystem(pmap,[1],[plane[2],0.0]) #take out the the reduced dimension (3d->2d)
grid_size = 2^5
#grid_edges::Vector{Float64} = [-9,0.0,-3.0,0.05]
#grid_edges::Vector{Float64} = [-10.0,-0.001,-2.0,0.5]
grid_edges::Vector{Float64} = [-10.0,-2.0] #grid on only 1d (x component)
Δt = 0.01
rw_ensemble = 1000 
nr_steps = 10000
writedlm(data_dir*"method_params.txt",[T,Ttr,grid_size,grid_edges])

orders::Vector{Int64} = [1,4,8,12] #[1,2,3,4] #[1,4,8,12]
const n::Int64 = 100 
const N_steps::Int64 = 500
const ensemble::Int64 = 100 #1000
const ϵ::Float64 = 1e-3 #1e-5

out_file_entropy = data_dir*"roessler_entropies"  * "_n_poinc_$T" * "Ttr_$Ttr" * "_plane_$plane_string" * "_b_0.3" * "_grid_$grid_size" * "_param_$b_start" * "_$b_stop" * ".txt"
out_file_lambda = data_dir*"roessler_lambdas" *  "_n_poinc_$T"* "Ttr_$Ttr" * "plane_$plane_string" * "_b_0.3" * "_grid_$grid_size" * "param_$b_start" * "_$b_stop" * ".txt"

function avg_crossing_time(pmap,T;Ttr)
	t_avg = 0.0

	#transient 
	for _ in 1:Ttr	
		step!(pmap)
	end

	t0 = current_crossing_time(pmap) #first current crossing time

	for _ in 1:T	
		step!(pmap)
    	t = current_crossing_time(pmap)
        t_avg += t - t0
        t0 = t
	end
	return t_avg/T
end


#---------------orbit diagram-----------------------

println("Orbit diagram calculation starting...")
od = orbitdiagram(pmap, 1, 2, b_vals; n=nr_points, Ttr=Ttr, show_progress=true)
writedlm(data_dir*"od_roessler_b_saved_z_T_$T"*"_Ttr_$Ttr"*"_b_$(b_start)_$b_stop"*".txt",stack(od)')
println("Done.")

#---------------normal lyapunovs--------------------
#calc here avg crossing time of pmap

println("Lyapunov exponent calculation starting...")
lyapunov_exponents = zeros(length(b_vals))
ds_copy = deepcopy(ds) #create independent ds 
pmap = PoincareMap(ds_copy,plane; rootkw=(xrtol=1e-10, atol=1e-10), direction=1) #becomes discrete system

avg_crossing_times = zeros(length(b_vals))
for (i,b) in enumerate(b_vals)
	@show b

	#avg crossing time
	set_parameter!(pmap,2,b)
	avg_crossing_times[i] = avg_crossing_time(pmap,1000;Ttr=1000)

	
	#lyap exp
	set_parameter!(ds_copy,2,b)
	lyap_exp = lyapunov(ds_copy,T;Ttr = Ttr)
	lyapunov_exponents[i] = lyap_exp
	
end

writedlm(data_dir*"lyapexps_roessler_b_T_$T"*"_Ttr_$Ttr"*"_b_$(b_start)_$b_stop"*".txt",hcat(b_vals,lyapunov_exponents))
writedlm(data_dir*"avg_cross_times_roessler_b_T_$T"*"_Ttr_$Ttr"*"_b_$(b_start)_$b_stop"*".txt",hcat(b_vals,avg_crossing_times))

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
	dim=1,
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

#writedlm("data/data_netmeasures_roessler_b_T_$T"*"_Ttr_$Ttr"*"_traj_ens"*"_rwens_$rw_ensemble"*"_nsteps_$nr_steps"*"_grid_$grid_size"*"standardpsos" *"_b_$(b_start)_$b_stop"*".txt",hcat(b_vals,analytic_entropies,analytic_lyapunovs,num_entropies,num_lyapunovs))



