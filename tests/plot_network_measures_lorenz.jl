using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using GraphPlot
using Random

# Parameters, preliminary data for S_KS and Lyapunov measure plots
grid = 20;
plane = (1,15.0);
rho_ensemble = 20; # Number of STNs for each value of ρ
rho = 180:0.003:182;
ensemble = 100; # Number of random walks
N_steps = 10000; # Length of random walks
L_measure = zeros(length(rho));
SKe = zeros(length(rho));
SK_analytic = zeros(length(rho));


@time for r in eachindex(rho)
    ent = 0.0;
    lya = 0.0;
    ent_analytic = 0.0;
    for i in 1:rho_ensemble
        @show i
        u0 = rand(Float64,3).*20 .-10
        system = Systems.lorenz(u0; ρ=rho[r]);
        psection = poincaresos(system, plane, 3000; Ttr=1000, direction=+1, rootkw = (xrtol = 1e-8, atol = 1e-8) );
        timeseries = psection[:,2:end]
        discrete_timeseries, vertex_names = timeseries_to_grid(timeseries, grid);
        stn_Q, stn_P = create_STN(discrete_timeseries, vertex_names)
        if is_strongly_connected(stn_P)
            result = walk_statistics(ensemble, stn_P, N_steps)
            ent += result[1]
            lya += result[2]
            ent_analytic += sinai_kolmogorov_entropy(stn_Q; graph_P=stn_P)
        end
    end
    SKe[r] = ent/rho_ensemble
    L_measure[r] = lya/rho_ensemble
    SK_analytic[r] = ent_analytic/rho_ensemble
end

# S_KS entropy plots
plot(rho, SKe,
ylims = (0,2),
yticks = 0:1:2,
xticks = 180:0.25:182,
linewidth = 3,
linecolor = :black,
xlabel = "ρ",
ylabel = "S_KS",
label = "Random path",
title = "20 initial conditions"
)

plot!(rho, SK_analytic,
ylims = (0,2),
yticks = 0:1:2,
xticks = 180:0.25:182,
linewidth = 1,
linecolor = :red,
xlabel = "ρ",
ylabel = "S_KS",
label = "Analyitic"
)

# Lyapunov measure plot
plot(rho, L_measure,
ylims = (0,3),
yticks = 0:1:3,
xticks = 180:0.25:182,
linewidth = 1,
linecolor = :black,
xlabel = "ρ",
ylabel = "Λ",
label = "Lyapunov measure",
title = "20 initial conditions"
)