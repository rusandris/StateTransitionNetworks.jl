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

Δt = 0.001;
plane = (1,15.0);
grid = 20;
ρ= 180.78;
u0 = rand(Float64,3).*50 .-25;
T = 5000;
ensemble = 10000
N_steps = 10000

system = Systems.lorenz(u0; ρ=ρ);
timeseries = trajectory(system, T; Δt=Δt, Ttr=500);

psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
d_traj, v_names = timeseries_to_grid(psection, grid);
stn, ret_code = create_stn(d_traj, v_names);
ret_code

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

L, S, Λ = walk_length_distribution(stn, ensemble, N_steps)

histogram(L, normalize=:pdf, xlabel="L", ylabel="p(L)", label=L"\rho=180.78")

S_hist = mean(L)/N_steps
S
Λ_hist = var(L, corrected=false)/N_steps
Λ

μ = S_hist*N_steps
σ = sqrt(Λ_hist*N_steps)
L_range = collect(range(μ-4σ,μ+4σ,1000))
plot!(L_range, gauss.(L_range,μ,σ), linewidth=4, label="Gauss")
μ+σ

