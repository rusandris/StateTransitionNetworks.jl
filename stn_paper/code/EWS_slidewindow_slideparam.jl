using TransitionsInTimeseries
using StateTransitionNetworks
using ComplexityMeasures
using StatsBase
using SparseArrays
using DynamicalSystems
using DelimitedFiles
cd(@__DIR__)
include("pipeline_functions.jl")

#------------------------Logistic stuff-----------------------

#ds = Systems.logistic()

function sliding_logistic(u,p,n)
	ϵ,n_min,n_max = p
	x = u[2]*u[1]*(1-u[1]) 
	r = n_min-1 <= n < n_max ? u[2] + ϵ : u[2] #slide 
	return SVector{2}(x,r)
end

r_min::Float64 = 3.8
r_max::Float64 = 3.84
r_crit::Float64 = 3.82842
#r_min::Float64 = 3.9
#r_max::Float64 = 3.91
#r_crit::Float64 = 3.9056 #critical value for the static (non-sliding) case
n = 100000
u0 = [0.3,r_min] #r starts at r_min

ϵ::Float64 = (r_max-r_min)/n #increment in r parameter
t_crit::Float64 = (r_crit - r_min)/ϵ #timepoint where it reaches r_crit
n_min::Int64 = 5000 #start of param sliding
n_max::Int64 = n_min + n #transient plus sliding time
p::Vector{Float64} = [ϵ,n_min,n_max]
ds = DeterministicIteratedMap(sliding_logistic,u0,p)

#rs = [3.8:Δp:3.84;]
grid_size = 20 #32
grid_edges = [0.0,1.0]
lag_logistic = 1
lags = [1:20;]
window_size = 1000#Int(n/20)
nr_ics = 3 #nr of initial conditions

#window parameters
indicator_window = (width = window_size, stride = 1)

times = [n_min:n_max;]
#gives the last time step of each window
t_windows = windowmap(last, times; indicator_window...)

#data containers to save
measures_container = []
ts_container = zeros(length(times))
params_container = zeros(length(times))

#study fluctuation of results
#run for many initial conditions
for i in 1:nr_ics

    orbit,_ = trajectory(ds,n_max-1,[rand(),r_min]) 
    #----------------------TransitionsInTimeseries-------------------
    #IDEA: use windowmapping from TransitionInTimeseries with indicators from StateTransitonNetworks

    #take only x component
    timeseries = orbit[:,1][n_min:end]
    if i == 1
        ts_container .= timeseries 
        params_container .= orbit[:,2][n_min:end]
    end

    #indicator time series
    measures = windowmap(timeseries -> early_warning_signals(timeseries,grid_size,grid_edges,lags), timeseries; indicator_window...)
    M = stack(measures)'
    @show size(M)
    push!(measures_container,M)
end    

#-------------------------------save data-------------------------------------
#TODO: mkdir data_dir if it doesn't exist
data_dir = "../data/supplimentary/sliding_param_measures/"


#leave off transient
ts_data = hcat(times,params_container,ts_container)
writedlm(data_dir*"logistic_map_timeseries_slideparam_tslength$n_max" * "n_tr$n_min" * "_r_min$r_min" * "_r_max$r_max" * ".txt",ts_data)
#measures
Ss = hcat(t_windows,[measures_container[i][:,1] for i in 1:nr_ics]...)
Λs = hcat(t_windows,[measures_container[i][:,2] for i in 1:nr_ics]...)
vars = hcat(t_windows,[measures_container[i][:,7] for i in 1:nr_ics]...)
acs = hcat(t_windows,[measures_container[i][:,8] for i in 1:nr_ics]...)

#write out time stuff as well
writedlm(data_dir*"logistic_map_time_param_slideparam" * "_r_min$r_min" * "_r_max$r_max" * ".txt",[n_min,n_max,t_crit])
#write out method parameters
writedlm(data_dir*"logistic_map_method_param_slideparam" * "_r_min$r_min" * "_r_max$r_max" * ".txt",[grid_size,window_size])
#write out method parameters
writedlm(data_dir*"logistic_map_system_param_slideparam" * "_r_min$r_min" * "_r_max$r_max" * ".txt",[r_min,r_max,r_crit])

#"$n_max" * "_window$window_size" * "_grid$grid_size" 
writedlm(data_dir*"logistic_map_entropies_slideparam.txt",Ss)
writedlm(data_dir*"logistic_map_lambdas_slideparam.txt",Λs)
writedlm(data_dir*"logistic_map_vars_slideparam.txt",vars)
writedlm(data_dir*"logistic_map_acs_slideparam.txt",acs)