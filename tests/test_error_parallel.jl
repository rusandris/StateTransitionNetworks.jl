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
import StateTransitionNetworks.sinai_kolmogorov_entropy

# PARALLEL ====================================================================================================================
println("Packages imported.")
flush(stdout)
Δt = 0.001;
plane = (1,15.0);
grid = 20;
rho = [180.1, 180.7, 180.78];
labels = ["chaos1", "chaos2", "ppchaos"];
ensemble = 10^4;

# Equidistant on log scale
T = [
    10, 17.783, 31.623, 56.234,
    100, 177.83, 316.23, 562.34,
    1000, 1778.3, 3162.3, 5623.4,
    10000, 17783, 31623, 56234, 100000
];


function calc_lyapunov_measure(ρ, total_time, u0)
    ret_code_stn = :Fail
    system = Systems.lorenz(u0; ρ=ρ);
    timeseries = trajectory(system, total_time; Δt=Δt, Ttr=1000);
    psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
    d_traj, v_names = timeseries_to_grid(psection, grid);
    stn, ret_code_stn = create_stn(d_traj, v_names; make_ergodic=true,verbose=false);
    if ret_code_stn == :Success
        P = prob_matrix(stn);
        l = lyapunov_measure(P)
        s = sinai_kolmogorov_entropy(P)
        if l[4] == :Success && s[2] == :Success
            return l[1], s[1], ret_code_stn
        else
            print("Fail")
            return -1, -1, l[4] && s[2]
        end
    else
        return -1, -1, ret_code_stn
    end    
end

function generate_initial_states(ensemble, ρ)
    init = SharedMatrix{Float64}(ensemble, 3)
    for i in 1:ensemble
        ret_code_stn = :Fail
        system = Systems.lorenz([0, 0, 0]; ρ=ρ);
        while ret_code_stn != :Success
            u0 = rand(Float64, 3).*200 .-100;
            timeseries = trajectory(system, 10, u0; Δt=Δt, Ttr=1000);
            psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
            d_traj, v_names = timeseries_to_grid(psection, grid);
            stn, ret_code_stn = create_stn(d_traj, v_names; make_ergodic=true,verbose=false);
            if ret_code_stn == :Success
                init[i, :] = u0
            end
        end
    end
    return init
end

for label in labels
    try
        rm("lorenz_error_*_$label.txt")
    catch
        continue
    end
end

L = SharedArray{Float64}(ensemble);
S = SharedArray{Float64}(ensemble);

for i in eachindex(rho)
    ρ=rho[i]
    println("ρ=",ρ)
    flush(stdout)
    label = labels[i]
    file = open("lorenz_error_lyapunov_$label.txt", "w")
    file2 = open("lorenz_error_entropy_$label.txt", "w")
    close(file)
    close(file2)
    println("Generating initial states")
    flush(stdout)
    states = generate_initial_states(ensemble, ρ)
    for time in T
        println(time)
        flush(stdout)
        @threads for j in 1:ensemble
            L[j], S[j] = calc_lyapunov_measure(ρ, time, states[j, :])
        end
        file = open("lorenz_error_lyapunov_$label.txt", "a")
        writedlm(file, [time mean(L[L .!= -1]) std(L[L .!= -1])])
        close(file)
        file2 = open("lorenz_error_entropy_$label.txt", "a")
        writedlm(file2, [time mean(S[S .!= -1]) std(S[S .!= -1])])
        close(file2)
    end
end

# PLOTTING RESULTS
# for i in eachindex(rho)
#     data =  readdlm("lorenz_error_$(labels[i]).txt", '\t', Float64, '\n')
#     T = data[:, 1] 
#     L_mean = data[:, 2]
#     L_std = data[:, 3]
#     error_plot = scatter(T, real(L_std./L_mean), xlabel = L"T", ylabel = L"σ_{Λ}/⟨Λ⟩",  label="", title="ρ=$(rho[i])", xaxis=:log, yaxis=:log)
#     mean_plot = scatter(T, L_mean, yerr=L_std./L_mean, xaxis=:log, xlabel=L"T", yaxis=L"Λ", title="ρ=$(rho[i])", label="")
#     savefig(error_plot, "lorenz_error_$(labels[i]).svg")
#     savefig(mean_plot, "lorenz_mean_$(labels[i]).svg")
# end
