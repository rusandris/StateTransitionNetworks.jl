using DelimitedFiles
using Printf
using DynamicalSystems
using StatsBase
#using ComplexityMeasures
cd(@__DIR__)
include("functions_chaotic_maps.jl")

#calc and save var and acf for Henon and logistic maps

const T::Int64 = 10^6 #10^8 
const T_string::String = @sprintf "%.E" T
const Ttr::Int64 = 10^3 #10^6 
const Ttr_string::String = @sprintf "%.E" Ttr

data_dir = "../data_test/supplimentary/static_measures_var_ac/"
mkpath(data_dir)

#TODO:speed this up!
#------------------------------Henon stuff--------------------------------

ds = Systems.henon()
Δp = 1e-3 #2*1e-4
as = [1.2:Δp:1.4;]
EWS_standard_henon = zeros(length(as),4)
orbits_henon = zeros(T,2) 
lag = 1

for (i,a) in enumerate(as)
    set_parameter!(ds,1,a)
    #@time traj,t = trajectory(ds,T;Ttr=Ttr)
    #@time traj = Matrix(traj)

    orbit!(orbits_henon,ds,T,Ttr;reinit=true)


    #autocorrelation and variance
    vs = vec(var(orbits_henon;dims=1,corrected=false)) #lag1, normalized by variance
    ac = vec(autocor(orbits_henon,[1])) #lag1, normalized by variance

    EWS_standard_henon[i,1:2] .= vs 
    EWS_standard_henon[i,3:4] .= ac 

end

writedlm(data_dir*"var_acf_henon_T$(T_string)_TTr$(Ttr_string)_rs_3.8_4.0.txt",EWS_standard_henon)

#------------------------Logistic stuff-----------------------

ds = Systems.logistic()
rs = [3.8:Δp:4.0;]
EWS_standard_log = zeros(length(as),2)
orbits_log = zeros(T,1) 
lag_logistic = 1

for (i,r) in enumerate(rs)
    set_parameter!(ds,1,r)
    #traj,t = trajectory(ds,T;Ttr=Ttr)
    #traj = Matrix(traj)

    orbit!(orbits_log,ds,T,Ttr;reinit=true)

    #autocorrelation and variance
    vs = var(orbits_log;corrected=false) #lag1, normalized by variance
    ac = autocor(orbits_log,[1])[1] #lag1, normalized by variance

    EWS_standard_log[i,1] = vs 
    EWS_standard_log[i,2] = ac 

end

writedlm(data_dir*"var_acf_log_T$(T_string)_TTr$(Ttr_string)_rs_3.8_4.txt",hcat(rs,EWS_standard_log))

#=
#check values on inset with finer param resolution

#------------------------------Henon stuff--------------------------------

ds = Systems.henon()
#as = [1.2:Δp:1.4;]
as = [1.305:Δp/4:1.32;]
EWS_standard_henon = zeros(length(as),4)
orbits_henon = zeros(T,2) 
lag = 1

for (i,a) in enumerate(as)
    @show a
    set_parameter!(ds,1,a)
    #@time traj,t = trajectory(ds,T;Ttr=Ttr)
    #@time traj = Matrix(traj)

    orbit!(orbits_henon,ds,T,Ttr;reinit=true)


    #autocorrelation and variance
    vs = vec(var(orbits_henon;dims=1,corrected=false)) #lag1, normalized by variance
    ac = vec(autocor(orbits_henon,[1])) #lag1, normalized by variance

    EWS_standard_henon[i,1:2] .= vs 
    EWS_standard_henon[i,3:4] .= ac 

end

writedlm("var_acf_henon_zoom_T1e8_TTr_1e6_as_$(as[1])_$(as[end])_4.txt",hcat(as,EWS_standard_henon))

#------------------------Logistic stuff-----------------------

ds = Systems.logistic()
#rs = [3.8:Δp:4.0;]
rs = [3.92:Δp/4:3.928;]
EWS_standard_log = zeros(length(rs),2)
orbits_log = zeros(T,1) 
lag_logistic = 1

for (i,r) in enumerate(rs)
    @show r
    set_parameter!(ds,1,r)
    #traj,t = trajectory(ds,T;Ttr=Ttr)
    #traj = Matrix(traj)

    @time orbit!(orbits_log,ds,T,Ttr;reinit=true)

    #autocorrelation and variance
    vs = var(orbits_log;corrected=false) #lag1, normalized by variance
    ac = autocor(orbits_log,[1])[1] #lag1, normalized by variance

    EWS_standard_log[i,1] = vs 
    EWS_standard_log[i,2] = ac 

end

writedlm("var_acf_log_zoom_T1e8_TTr_1e6_rs_$(rs[1])_$(rs[end]).txt",hcat(rs,EWS_standard_log))

=#