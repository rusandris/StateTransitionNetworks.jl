using StateTransitionNetworks
using DynamicalSystems
using ChaosTools
using Plots
using LinearAlgebra
using StatsBase
using LaTeXStrings
using DelimitedFiles

import StateTransitionNetworks.lyapunov_measure

Δt = 0.001;
plane = (1,15.0);
grid = 20;
rho = [180.1, 180.7, 180.78];
labels = ["chaos1", "chaos2", "ppchaos"];
ensemble = 10^5;

# Equidistant on log scale
T = [
    10, 17.783, 31.623, 56.234,
    100, 177.83, 316.23, 562.34,
    1000, 1778.3, 3162.3, 5623.4,
    10000, 17783, 31623, 56234, 100000
];
for i in eachindex(rho)
    @show ρ=rho[i]
    label = labels[i]
    L_data = [];
    L_value = [];
    for t in eachindex(T)
        @show T[t]
        L = []
        local P, Q
        while length(L) < ensemble
            u0 = rand(Float64,3).*50 .-25;
            system = Systems.lorenz(u0; ρ=ρ);
            timeseries = trajectory(system, T[t]; Δt=Δt, Ttr=1000);
            psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
            d_traj, v_names = timeseries_to_grid(psection, grid);
            stn, ret_code_stn = create_stn(d_traj, v_names, make_ergodic=true);
            P = prob_matrix(stn);
            Q = weight_matrix(stn);
            if !any(isnan, P)
                lyapunov = lyapunov_measure(P)[1]
                push!(L, lyapunov)
            end
        end
        push!(L_data, std(L))
        push!(L_value, mean(L))
    end        
    file = open("lorenz_error_$label.txt", "w")
    writedlm(file, [T real(L_value) real(L_data)])
    close(file)
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