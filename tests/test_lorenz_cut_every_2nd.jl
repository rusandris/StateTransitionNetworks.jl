using StateTransitionNetworks
using DynamicalSystems
using ChaosTools
using Plots
using StatsPlots
import StateTransitionNetworks.lyapunov_measure
import StateTransitionNetworks.sinai_kolmogorov_entropy
using LaTeXStrings
using DelimitedFiles

# ===============================================================================================
# CUTTING TIMESERIES INTO SEGMENTS, DISCARDING EVERY 2ND SEGMENT
# NUMBER OF CUTS IS INCREASED AND NETWORK MEASURES ARE CALCULATED
# ===============================================================================================

Δt = 0.001;
plane = (1,15.0);
grid = 20;
rho = [180.1, 180.7, 180.78, 181.1];
labels = ["chaos1", "chaos2", "ppchaos", "periodic"];
T = 13107.2;

cuts = [2, 2^2, 2^4, 2^6, 2^8, 2^10, 2^12, 2^14];

n_ens = 100;

for index in eachindex(rho)
    label = labels[index]
    @show label
    ρ = rho[index]
    @time for j in 1:n_ens
        @show j
        stream = open("lorenz_cut_$label.txt", "a")
        u0 = rand(Float64,3).*50 .-25;
        system = Systems.lorenz(u0; ρ=ρ);
        timeseries,  = trajectory(system, T; Δt=Δt, Ttr=500)[1:end-1];
        L = length(timeseries);
        timeseries = Matrix(timeseries);
        S_data = []
        L_data = []
        V_data = []
        C_data = []
        for n in eachindex(cuts)
            len = div(L,cuts[n]+1)
            new_timeseries = timeseries[1:len, :]
            for i in 1:cuts[n]
                if i % 2 == 0
                    new_timeseries = vcat(new_timeseries, timeseries[(len*i+1):(len*(i+1)), :])
                end
            end
            new_timeseries = Dataset(new_timeseries)
            psection = DynamicalSystemsBase.poincaresos(new_timeseries, plane; direction=+1, save_idxs=[2,3])
            d_traj, v_names = timeseries_to_grid(psection, grid);
            stn, ret_code = create_stn(d_traj, v_names; make_ergodic=true, verbose=true);
            if ret_code == :Success
                P = prob_matrix(stn);
                push!(S_data, sinai_kolmogorov_entropy(P))
                l, v, c, success = lyapunov_measure(P)
                push!(L_data, l)
                push!(V_data, v)
                push!(C_data, 2*c)
            end
        end
        writedlm(stream, [cuts S_data L_data V_data C_data])
        close(stream)
    end
end
# ===========================================================================================
# CALCULATION ENDS HERE
# EVERYTHING IS WRITTEN OUT TO FILES
# ===========================================================================================


# ===========================================================================================
# PLOTTING STARTS HERE
# ===========================================================================================

Δt = 0.001;
plane = (1,15.0);
grid = 20;
rho = [180.1, 180.7, 180.78, 181.1];
labels = ["chaos1", "chaos2", "ppchaos", "periodic"];
T = 13107.2;
max_period = 1.024;

# CALCULATING REFERENCE VALUES
measures = zeros(Float64, (4,4));
ensemble = 50;
for index in eachindex(rho)
    ρ = rho[index]
    @show ρ
    entropy = 0.0
    lyapunov = 0.0
    variance = 0.0
    covariance = 0.0
    for i in 1:ensemble
        u0 = rand(Float64,3).*50 .-25;
        system = Systems.lorenz(u0; ρ=ρ);
        timeseries,  = trajectory(system, T; Δt=Δt, Ttr=500)[1:end-1];
        psection = DynamicalSystemsBase.poincaresos(timeseries, plane; direction=+1, save_idxs=[2,3])
        d_traj, v_names = timeseries_to_grid(psection, grid);
        stn, ret_code = create_stn(d_traj, v_names);
        if ret_code == :Success
            P = prob_matrix(stn)
            entropy += sinai_kolmogorov_entropy(P)
            l, v, c, success = lyapunov_measure(P)
            lyapunov += l
            variance += v
            covariance += 2*c
        else
            i = i-1
        end
    end
    measures[index, :] = [entropy, lyapunov, variance, covariance]./ensemble
end

measures

# ACTUAL PLOTTING

for index in eachindex(rho)
    label = labels[index]
    ρ = rho[index]
    data = readdlm("lorenz_cut_$label.txt", '\t', Float64, '\n')

    cuts = data[:, 1]
    sort!(cuts)
    unique!(cuts)
    S = data[:, 2]
    L = data[:, 3]
    V = data[:, 4]
    C = data[:, 5]
    entropy, lyapunov, variance, covariance = measures[index, :]

    ticks_place = log.(2, T ./(cuts .+ 1))
    ticks_place[1]=13
    ticks = T ./(cuts .+ 1)
    ticks = round.(ticks, digits=1)

    splot = hline([entropy], linewidth = 3, linecolor=:red, label="Reference value")
    boxplot!( splot,
        ticks_place, S, 
        xlabel= L"\Delta T",
        ylabel = L"S_{KS}",
        xticks = (ticks_place, ticks),
        title = "ρ=$ρ, T=$T",
        label = "",
        fillalpha=0.75,
        fillcolor=:skyblue2)

    vline!(splot,[max_period], linewidth = 3, linecolor=:green, label=L"\Delta t_P^{(max)}=1.024", legend=:best)
    savefig(splot, "lorenz_cut_plot_entropy_$label.svg")

    lplot = hline([lyapunov], linewidth = 3, linecolor=:red, label="Reference value")
    boxplot!( lplot,
        ticks_place, L, #L_matrix[9:end], 
        xlabel= L"\Delta T",
        xticks = (ticks_place, ticks),
        ylabel = L"Λ",
        title = "ρ=$ρ, T=$T",
        label = "",
        fillalpha=0.75,
        fillcolor=:gray69)
    vline!(lplot,[max_period], linewidth = 3, linecolor=:green, label=L"\Delta t_P^{(max)}=1.024", legend=:best)
    savefig(lplot, "lorenz_cut_plot_lyapunov_$label.svg")

    lyap = boxplot(
        ticks_place, L, 
        xlabel= L"\Delta T",
        xticks = (ticks_place, ticks),
        ylabel = "",
        title = "ρ=$ρ, T=$T",
        label = L"Λ",
        fillalpha=0.75,
        fillcolor=:gray69)
    boxplot!( lyap,
        ticks_place, V, 
        label = L"\sigma^2",
        fillalpha=0.75,
        fillcolor=:tomato)
    boxplot!( lyap,
        ticks_place, C, 
        label = L"2Cov",
        fillalpha=0.75,
        fillcolor=:palegreen)
    hline!(lyap, [lyapunov], linewidth = 3, linecolor=:red, label="Reference value")
    vline!(lyap,[max_period], linewidth = 3, linecolor=:green, label=L"\Delta t_P^{(max)}=1.024", legend=:best)
    savefig(lyap, "lorenz_cut_plot_lyapunov_co_variance_$label.svg")
end