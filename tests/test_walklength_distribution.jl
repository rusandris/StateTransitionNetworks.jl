using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using GraphPlot
using Random
using Statistics
import StatsBase: sample, mean, var, Weights
using LaTeXStrings
using LinearAlgebra

function walk_length_distribution(stn, ensemble, N_steps)
    walk_length = Vector{Float64}(undef, ensemble)
    for i in 1:ensemble
        walk_length[i] = random_walk_on_stn(stn, N_steps)
    end
   	entropy = mean(walk_length)/N_steps
    lyapunov_measure = var(walk_length,corrected=false)/N_steps
    return walk_length, entropy, lyapunov_measure
end

gauss(x,μ,σ) = 1/(σ*sqrt(2π)) * exp(-(x - μ)^2/(2σ^2)) 


Δt = 0.001;
plane = (1,15.0);
grid = 20;
u0 = rand(Float64,3).*50 .-25;
T = 5000;
ensemble = 10000
N_steps = 10000

rho_arr = [180.1, 180.7, 180.78]
rho_arr = [180.1]
color_arr = [:orange, :red, :green]
pl = plot()
for (i,ρ) in enumerate(rho_arr)
    system = Systems.lorenz(u0; ρ=ρ);
    timeseries = trajectory(system, T; Δt=Δt, Ttr=500);
    psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
    d_traj, v_names = timeseries_to_grid(psection, grid);
    stn, ret_code = create_stn(d_traj, v_names);
    ret_code
    L, S, Λ = walk_length_distribution(stn, ensemble, N_steps)
    histogram!(pl, L, normalize=:pdf, bins=100,
        color=color_arr[i],
        alpha=0.3,
        lw = 0,
        xlabel=L"L",
        ylabel=L"$p(L)$",
        xguidefontsize=18,
        yguidefontsize=18,
        tickfontsize=10,
        #label=L"\rho_{%$(i)}=%$(ρ)",
        label=L"\rho=%$(ρ)",
        legendfontsize=15,
        legend=:best,
        dpi=300
        )

    S_hist = mean(L)/N_steps
    Λ_hist = var(L, corrected=false)/N_steps

    μ = S_hist*N_steps
    σ = sqrt(Λ_hist*N_steps)
    L_range = collect(range(μ-4σ,μ+4σ,1000))
    plot!(pl, L_range, gauss.(L_range,μ,σ),
        color=color_arr[i],
        linewidth=2,
        xlabel=L"L",
        ylabel=L"$p(L)$",
        xguidefontsize=18,
        yguidefontsize=18,
        tickfontsize=10,
        label=L"\mathcal{N}(S\,t,\Lambda\,t)",
        #label=nothing,
        legendfontsize=15,
        legend=:best,
        dpi=300
        ) 
    savefig(pl, "./tests/walk_length_distribution_rho=$(ρ)_Nsteps=$(N_steps)_ensemble=$(ensemble).svg")

end
plot!(pl, xlim=(4000,12500), size = (700,380), pad=10)
savefig(pl, "./tests/walk_length_distribution_Nsteps=$(N_steps)_ensemble=$(ensemble).svg")




μ+σ

