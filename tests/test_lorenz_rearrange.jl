using Random
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
# CUTS WITH REORDERING
# ===============================================================================================


Δt = 0.001;
plane = (1,15.0);
grid = 20;
rho = [180.1, 180.7, 180.78, 181.1];
labels = ["chaos1", "chaos2", "ppchaos", "periodic"];
T = 13107.2;
max_period = 1.024;
cuts = [2, 2^2, 2^4, 2^6] #, 2^8, 2^10, 2^12, 2^14];
n_ens = 5;


for index in eachindex(rho)
    label = labels[index]
    @show label
    ρ = rho[index]
    for j in 1:n_ens
        @show j
        u0 = rand(Float64,3).*50 .-25;
        system = Systems.lorenz(u0; ρ=ρ);
        timeseries = trajectory(system, T; Δt=Δt, Ttr=500)[1:end-1];
        L = length(timeseries);
        stream = open("lorenz_rearrange_$label.txt", "a")
        S_data = []
        L_data = []
        V_data = []
        C_data = []
        for n in eachindex(cuts)
            len = div(L,cuts[n]+1)

            u = timeseries[len,:]
            system = Systems.lorenz(u; ρ=ρ);
            new_timeseries = Matrix(trajectory(system, len*Δt; Δt=Δt))
            for i in shuffle(2:cuts[n])
                u = timeseries[i*len,:]
                system = Systems.lorenz(u; ρ=ρ);
                new_timeseries = vcat(new_timeseries, Matrix(trajectory(system, len*Δt; Δt=Δt)))
            end

            new_timeseries = Dataset(new_timeseries)
            psection = ChaosTools.poincaresos(new_timeseries, plane; direction=+1, idxs=[2,3])
            d_traj, v_names = timeseries_to_grid(psection, grid);
            stn, ret_code = create_stn(d_traj, v_names; make_ergodic=true, verbose=true)
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
        @show i
        u0 = rand(Float64,3).*50 .-25;
        system = Systems.lorenz(u0; ρ=ρ);
        timeseries = trajectory(system, T; Δt=Δt, Ttr=500)[1:end-1];
        psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3])
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
    measures[1, :] = [entropy, lyapunov, variance, covariance]./ensemble
end

measures



# ACTUAL PLOTTING

for index in eachindex(rho)
    label = labels[index]
    ρ = rho[index]
    data = readdlm("lorenz_rearrange_$label.txt", '\t', Float64, '\n')

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
    savefig(splot, "lorenz_rearrange_plot_entropy_$label.svg")

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
    savefig(lplot, "lorenz_rearrange_plot_lyapunov_$label.svg")

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
    savefig(lyap, "lorenz_rearrange_plot_lyapunov_co_variance_$label.svg")
end