using Revise
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

import StateTransitionNetworks.lyapunov_measure
import SparseArrays: spzeros 
import StateTransitionNetworks.renormalize
import Graphs: DiGraph

function timeseries_to_common_grid(timeseries, grid, x_min, x_max, y_min, y_max)    
    T = length(timeseries[:,1])
    dx = 0.5*(x_max-x_min)/grid
    dy = 0.5*(y_max-y_min)/grid
    x_grid = range(x_min-dx, x_max+dx, grid);
    x_min = x_grid[1]
    y_grid = range(y_min-dy, y_max+dy, grid);
    y_min = y_grid[1]
    x_n = Vector{Int64}(undef, T)
    y_n = Vector{Int64}(undef, T)
    x_n, y_n = [], [];

    for row in eachrow(timeseries)
        y = floor(Int,(row[2]-y_min)/Float64(y_grid.step))+1
        x = floor(Int,(row[1]-x_min)/Float64(x_grid.step))+1
        push!(x_n, x)
        push!(y_n, y)
    end
    d_timeseries = [x_n y_n]
    return d_timeseries
end



function add_timeseries(dt_list, grid;make_ergodic=make_ergodic,verbose=verbose)
    Q_null = zeros(Int32, grid*grid, grid*grid) # Null matrix with all possible transitions
    vertex_names = [];  # Future name of vertices
    vertex_place = [];  # Rows, and column number in Q_null for a given vertex
    M = zeros(grid,grid)
    nr_vertices = 0;
    for dt in dt_list
        states = []
        O = zeros(Int32, grid*grid, grid*grid)
        for row in eachrow(dt)
            v = (row[1]-1)*grid + row[2]
            x, y = row
            if M[y,x] == 0
                nr_vertices += 1 
                push!(vertex_names, [nr_vertices, x, y])
                push!(vertex_place, v)
                M[y,x] = 1
            end    
            push!(states, v) 
        end
        #count transitions
        next_states = circshift(states,-1)
        for i in eachindex(states[1:end-1])
            O[states[i],next_states[i]] += 1
        end
        Q_null = Q_null + O
    end
    vertex_names = reduce(hcat, vertex_names)'
    Q = Q_null[vertex_place, vertex_place]
	#weight and transition probability matrices
    P = spzeros(nr_vertices, nr_vertices)
	#normalize Q and fill P by normalizing rows
    Q = Q./sum(Q)

	P = renormalize(Q)
	#create directed metagraph with static label and metadata types and default weight 0
    stn, ret_code = create_stn(P;make_ergodic=make_ergodic,verbose=verbose)
    for v in 1:nr_vertices
        x,y = vertex_names[v,2:end] 
        stn[v] = Dict(:x => x,:y => y)
    end
    return stn, ret_code
end

function get_grid_edges(psections)
    x_max, y_max = maximum(reduce(vcat, maximum.(Matrix.(psections); dims=1)); dims=1)
    x_min, y_min = minimum(reduce(vcat, minimum.(Matrix.(psections); dims=1)); dims=1)
    return x_min, x_max, y_min, y_max
end

# ===============================================================================================
# CUTS WITH REORDERING
# ===============================================================================================


Δt = 0.001;
plane = (1,15.0);
grid = 20;
rho = [180.1, 180.7, 180.78];
labels = ["chaos1", "chaos2", "ppchaos"];
T = 13107.2;
max_period = 1.024;
cuts = [2, 2^2, 2^4, 2^6, 2^8, 2^10, 2^12, 2^14];
n_ens = 10;

for label in labels
    try
        rm("lorenz_rearrange_$label.txt")
    catch
        continue
    end
end


@time for index in eachindex(rho)
    label = labels[index]
    @show label
    ρ = rho[index]
    for j in 1:n_ens
        @show j
        u0 = rand(Float64,3).*50 .-25;
        system = Systems.lorenz(u0; ρ=ρ);
        timeseries = trajectory(system, T; Δt=Δt, Ttr=500)[1:end-1];
        psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
        L = length(timeseries);
        stream = open("lorenz_rearrange_$label.txt", "a")
        S_data = []
        L_data = []
        V_data = []
        C_data = []
        for n in eachindex(cuts)
            len = div(L,cuts[n]+1)
            psection_list = []
            for i in shuffle(1:cuts[n])
                u = timeseries[Int(i*len),:]
                system = Systems.lorenz(u; ρ=ρ);
                new_timeseries = trajectory(system, len*Δt; Δt=Δt)
                psection = ChaosTools.poincaresos(new_timeseries, plane; direction=+1, idxs=[2,3]);
                if length(psection) != 0
                    push!(psection_list, psection)
                end
            end
            x_min, x_max, y_min, y_max = get_grid_edges(psection_list)
            dt_list = []
            for psection in psection_list
                dt = timeseries_to_common_grid(psection, grid, x_min, x_max, y_min, y_max);
                push!(dt_list, dt)
            end

            stn, ret_code = add_timeseries(dt_list, grid; make_ergodic=true, verbose=false)
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
    measures[index, :] = [entropy, lyapunov, variance, covariance]./ensemble
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
    vline!(lyap,[max_period], linewidth = 3, linecolor=:green, label=L"\Delta t_P^{(max)}=1.024", legend=:bottomright)
    savefig(lyap, "lorenz_rearrange_plot_lyapunov_co_variance_$label.svg")
end