using StateTransitionNetworks
using DynamicalSystems
using ChaosTools
using DelimitedFiles
using SparseArrays
using Printf
include("functions_chaotic_maps.jl")

println("Packages loaded...")
println("Number of threads available: ",Threads.nthreads())

flush(stdout)

#--------------------parameters--------------------

#for broader control parameter range
Δp::Float64 = 1e-3 	#2*1e-4
ps::Vector{Float64} = [3.5:Δp:4.0;]

const T::Int64 = 10^6 #10^8
const T_string::String = @sprintf "%.E" T
const Ttr::Int64 = 10^3 #10^6 
const Ttr_string::String = @sprintf "%.E" Ttr
const grid_size::Int64 = 10 #32
const grid_edges::Vector{Float64} = [0.0,1.0]

const orders::Vector{Int64} = [1,2,3,4] #[1,4,8,12]
const n::Int64 = 100
const N_steps::Int64 = 500
const ensemble::Int64 = 100 #1000
const ϵ::Float64 = 1e-3 #1e-5
const u0::Vector{Float64} = [0.3]

const ps_renyi::Vector{Float64} = [3.82842,3.835,3.85,4.0]
const qs::Vector{Float64} = [0.0:0.01:2;]


data_dir = "../data/supplimentary/logistic_data_sm/" 
mkpath(data_dir)

out_file_entropy = data_dir*"logistic_entropies"  * "_T$T_string" * "_Ttr$Ttr_string" * "_grid_$grid_size" * "param_$(ps[1])" * "_$(ps[end])" * ".txt"
out_file_lambda = data_dir*"logistic_lambda"  * "_T$T_string" * "_Ttr$Ttr_string" * "_grid_$grid_size" * "param_$(ps[1])" * "_$(ps[end])" * ".txt"


global ds = Systems.logistic()

#--------------------------------------renyi entropy spectrums------------------------------------
writedlm(data_dir * "renyi_param_values" * ".txt",ps_renyi)
calc_save_renyi_spectrums_grid(ds,qs,ps_renyi,orders;param_index=1,grid_size=grid_size,u0=[0.3],grid_edges=grid_edges,reinit=true,T=T,Ttr=Ttr,out_dir=data_dir)

@info "Renyi entropy spectrum done."
flush(stderr)

#------------------------------------------orbit diagram-----------------------------------------

od = orbitdiagram(ds,1,1,ps; n = n,Ttr=Ttr)
writedlm(data_dir * "logistic_od" * "_T$T_string" * "_Ttr$Ttr_string" * "_grid_$grid_size" * "_nr_param_$(length(ps))" * "param_$(ps[1])" * "_$(ps[end])" * ".txt",od)

@info "Orbit diagram done."
flush(stderr)

#----------------------------------------------calc_measures-----------------------------------------

writedlm(data_dir * "orders" * ".txt",orders)
writedlm(data_dir * "logistic_p_values" * "_nr_param_$(length(ps))"*".txt",ps)

@info "Calculating measures...."
flush(stderr)

calc_measures(ds,ps,orders;
    out_file_entropy=out_file_entropy,
    out_file_lambda=out_file_lambda,
    sts_eltype=Int128,
    param_index=1,
    grid_size=grid_size,
    grid_edges=grid_edges,
    u0=u0,
    T=T,
    Ttr=Ttr,
    reinit=true,
    ensemble = ensemble,
    N_steps = N_steps,
    ϵ=ϵ)

@info "Measure calculation done and saved."
flush(stderr)

#------------------------------------calc lyapunov exponents (LLE)-----------------------------------

lyap_exps = calc_lyapunovs(ds,ps,5000,5000;u0=u0)
writedlm(data_dir*"lyapexps" * "_T$T_string" * "_Ttr$Ttr_string" * "_grid_$grid_size" * "_nr_param_$(length(ps))" * "param_$(ps[1])" * "_$(ps[end])" * ".txt",hcat(ps,lyap_exps))

#-----------------------------------------------direct random walk simulation stats-----------------------------------------


const ps_log_rw::Vector{Float64} = [3.82842,3.85,4.0]
const order_rw::Int64 = 1

ds = Systems.logistic()

logistic_out1 = data_dir * "logistic_walk_lengths_ens_$(ensemble)" * "_t$(N_steps)" * "_T$T_string" * "_Ttr$Ttr_string" * "_order_$order_rw" * "_grid$(grid_size)" * ".txt"
logistic_out2 = data_dir * "logistic_wl_stats_ens_$(ensemble)" * "_t$(N_steps)" * "_T$T_string" * "_Ttr$Ttr_string" * "_order_$order_rw" * "_grid$(grid_size)" * ".txt"

walk_lengths_log,means_log,variances_log,Ss_log,Λs_log = walk_length_statistics(ds, ps_log_rw; grid_size=grid_size,order=order_rw, T = T, Ttr = Ttr, u0=[0.3],ensemble = ensemble, N_steps = N_steps)
writedlm(logistic_out1,walk_lengths_log)
writedlm(logistic_out2,hcat(means_log,variances_log,Ss_log,Λs_log))