using StateTransitionNetworks
using Plots
using DynamicalSystems
using SparseArrays
using StatsBase

using LinearAlgebra
using DelimitedFiles
using LaTeXStrings
using Printf
cd(@__DIR__)

figs_dir="../figs/"
mkpath(figs_dir)

data_dir="../data/supplimentary/crit_tent_maps/"
mkpath(data_dir)

include("functions_chaotic_maps.jl")


#-----------------------------------critical map--------------------------------------------

function f_crit(u, p, n)
    r = p[1]
    x = u[1]
    return SVector(1 - abs(x^r-(1-x)^r)^(1.0/r))
end


function static_renyi_entropy(ρ, q)
    Hq = 0.
    if isapprox(q,1.0)
        for i in 1:length(ρ)
            Hq += -ρ[i]*log(ρ[i])
        end
        return Hq
    else
        for i in 1:length(ρ)
            Hq += (ρ[i])^q
        end
        Hq = log(Hq)/(1-q)
        return Hq
    end
end

function static_renyi_entropy_spectrum(ρ, qs=collect(0:0.01:2))
    Hqs = zeros(length(qs))
    for (i,q) in enumerate(qs)
        Hqs[i] = static_renyi_entropy(ρ, q)
    end
    return Hqs
end

###############################
### Measures for a single value
###############################

r = 2.
ds_crit = DeterministicIteratedMap(f_crit, [0.4], [r])


T::Int64 = 10^7
T_string::String = @sprintf "%.E" T
Ttr::Int64 = 10^5
Ttr_string::String = @sprintf "%.E" Ttr
grid_edges::Vector{Float64} = [0,1]
writedlm(data_dir*"method_params.txt",[T,Ttr,grid_size])

timeseries, = trajectory(ds_crit, T, [0.3]; Ttr=Ttr);

# compute DYNAMICAL Renyi entropy spectrum for HIGHER ORDER states
grid_size = 2^10
qs = collect(0.01:0.01:2)
Os = [1,2,4,8,10]
writedlm(data_dir*"orders_crit.txt",Os)

sts = zeros(Int128,length(timeseries))
timeseries_to_grid!(sts, timeseries, grid_size; grid_edges = [0., 1.]);
sts2 = zeros(Int128, length(sts))
for (i,o) in enumerate(Os[1:end])
        sts2 .= sts
        @show o
        @time higher_order_symbolics!(sts2, o)
        nr_symbols = length(unique(sts2))
        @show nr_symbols
        @time P,Q,x = calculate_transition_matrix(@view sts2[1:end-o+1]);
        Hs = renyi_entropy_spectrum(P, qs; x=x, verbose=true)
        f_name = data_dir*"critical_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(Int(o)).dat"
        writedlm(f_name,hcat(qs, Hs))
end


#######################################
# compute STATIC Renyi entropy spectrum
#######################################

sts = timeseries_to_grid(timeseries, grid_size; grid_edges = [0., 1.]);
@time P,Q,x = calculate_transition_matrix(sts);
Hqs = static_renyi_entropy_spectrum(x, qs)
plot!(qs,Hqs)

#Os = [8,12,14,16,30]
Os = [1,2,4,8,10,12,14,16,30]

writedlm(data_dir*"orders_crit_static.txt",Os)

# higher order states
sts = zeros(Int128,length(timeseries))
timeseries_to_grid!(sts, timeseries, grid_size; grid_edges = [0., 1.]);
sts2 = zeros(Int128, length(sts))
for (i,o) in enumerate(Os[1:end])
        sts2 .= sts
        @show o
        @time higher_order_symbolics!(sts2, o)
        nr_symbols = length(unique(sts2))
        @show nr_symbols

        #use simple countmap to get state distribution
        #@time P,Q,x = calculate_transition_matrix(@view sts2[1:end-o+1]);
        @time x = collect(values(countmap(@view sts2[1:end-o+1]))) ./ length(sts2)

        Hs = static_renyi_entropy_spectrum(x, qs)
        f_name = data_dir*"critical_static_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(Int(o)).dat"
        writedlm(f_name,hcat(qs, Hs))
end

#-----------------------------------tent map--------------------------------------------

function f_tent(u, p, n)
    r = p[1]
    x = u[1]
    if x<=r
        return SVector(x/r)
    else
        return SVector((1-x)/(1-r))
    end
end

function analytical_renyi_entropy(w, q)
    if q==0
        return log(2)
    elseif q==1
        return -w*log(w)-(1-w)*log(1-w)
    else
        return log(w^q+(1-w)^q)/(1-q)
    end
end

function analytical_renyi_entropy_spectrum(w, qs; verbose=true)
    Hs = zeros(length(qs))
    for (i,q) in enumerate(qs)
        verbose && @show q
        Hs[i] = analytical_renyi_entropy(w, q)
    end
    return Hs
end

#define tent map
r = 0.8
ds_tent = DeterministicIteratedMap(f_tent, [0.4], [r])

grid_size = 2^5
qs = collect(0.01:0.01:2)
Os = [1,2,4,8,10]

timeseries, = trajectory(ds_tent, T, [0.3]; Ttr=Ttr);

sts = timeseries_to_grid(timeseries, grid_size; grid_edges = [0., 1.]);

# compute Renyi entropy spectrum for HIGHER ORDER states

writedlm(data_dir*"orders_tent_dynamic.txt",Os)
sts = zeros(Int128,length(timeseries))
timeseries_to_grid!(sts, timeseries, grid_size; grid_edges = [0., 1.]);
sts2 = zeros(Int128, length(sts))
for (i,o) in enumerate(Os[1:end])
        sts2 .= sts
        @show o
        @time higher_order_symbolics!(sts2, o)
        nr_symbols = length(unique(sts2))
        @show nr_symbols
        @time P,Q,x = calculate_transition_matrix(@view sts2[1:end-o+1]);
        Hs = renyi_entropy_spectrum(P, qs; x=x, verbose=true)
        f_name = data_dir*"asymmetric_triangular_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(Int(o)).dat"
        writedlm(f_name,hcat(qs, Hs))
end


#######################################
# compute STATIC Renyi entropy spectrum
#######################################
Os_static = [8,12,14,16,30]
writedlm(data_dir*"orders_tent_static.txt",Os_static)
sts = timeseries_to_grid(timeseries, grid_size; grid_edges = [0., 1.]);

# compute for higher order states
sts = zeros(Int128,length(timeseries))
timeseries_to_grid!(sts, timeseries, grid_size; grid_edges = [0., 1.]);
sts2 = zeros(Int128, length(sts))
for (i,o) in enumerate(Os_static[1:end])
        sts2 .= sts
        @show o
        @time higher_order_symbolics!(sts2, o)
        nr_symbols = length(unique(sts2))
        @show nr_symbols
        #@time P,Q,x = calculate_transition_matrix(@view sts2[1:end-o+1]);

        #use simple countmap to get state distribution
        #@time P,Q,x = calculate_transition_matrix(@view sts2[1:end-o+1]);
        @time x = collect(values(countmap(@view sts2[1:end-o+1]))) ./ length(sts2)

        Hs = static_renyi_entropy_spectrum(x, qs)
        f_name = data_dir*"asymmetric_triangular_static_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(Int(o)).dat"
        writedlm(f_name,hcat(qs, Hs))
end

