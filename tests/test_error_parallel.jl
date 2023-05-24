using StateTransitionNetworks
using DynamicalSystems
using ChaosTools
using Plots
using LinearAlgebra
using StatsBase
using LaTeXStrings
using DelimitedFiles
using SharedArrays
using .Threads


import StateTransitionNetworks.lyapunov_measure

# PARALLEL ====================================================================================================================

Δt = 0.001;
plane = (1,15.0);
grid = 20;
rho = [180.1, 180.7, 180.78];
labels = ["chaos1", "chaos2", "ppchaos"];
ensemble = 10^3;

# Equidistant on log scale
T = [
    10, 17.783, 31.623, 56.234,
    100, 177.83, 316.23, 562.34,
    #1000, 1778.3, 3162.3, 5623.4,
    #10000, 17783, 31623, 56234, 100000
];


function calc_lyapunov_measure(ρ, T)
    ret_code_stn = :Fail
    local stn, L
    while ret_code_stn != :Success
        u0 = rand(Float64,3).*50 .-25;
        system = Systems.lorenz(u0; ρ=ρ);
        timeseries = trajectory(system, T; Δt=Δt, Ttr=1000);
        psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
        d_traj, v_names = timeseries_to_grid(psection, grid);
        stn, ret_code_stn = create_stn(d_traj, v_names; make_ergodic=true,verbose=false);
        if ret_code_stn == :Success
            P = prob_matrix(stn);
            l = lyapunov_measure(P)
            if l[4] == :Success
                return l[1]
            else
                ret_code_stn = :Fail
            end
        end
    end    
end


for label in labels
    try
        rm("lorenz_error_$label.txt")
    catch
        continue
    end
end

L = SharedArray{Float64}(ensemble);
for i in eachindex(rho)
    @show ρ=rho[i]
    label = labels[i]
    L_mean = []
    L_std = []
    file = open("lorenz_error_$label.txt", "w")
    close(file)
    for time in T
        @threads for j in 1:ensemble
            L[j] = calc_lyapunov_measure(ρ, time)
            end
        push!(L_mean, mean(L))
        push!(L_std, std(L))
        file = open("lorenz_error_$label.txt", "a")
        writedlm(file, [time mean(L) std(L)])
        close(file)    
    end
end

# PLOTTING RESULTS
for i in eachindex(rho)
    data =  readdlm("lorenz_error_$(labels[i]).txt", '\t', Float64, '\n')
    T = data[:, 1] 
    L_mean = data[:, 2]
    L_std = data[:, 3]
    error_plot = scatter(T, real(L_std./L_mean), xlabel = L"T", ylabel = L"σ_{Λ}/⟨Λ⟩", xaxis=:log, label="", title="ρ=$(rho[i])")
    mean_plot = scatter(T, L_mean, yerr=L_std./L_mean, xaxis=:log, xlabel=L"T", yaxis=L"Λ", title="ρ=$(rho[i])", label="")
    savefig(error_plot, "lorenz_error_$(labels[i]).svg")
    savefig(mean_plot, "lorenz_mean_$(labels[i]).svg")
end