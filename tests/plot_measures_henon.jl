using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using GraphPlot
using Random
using Statistics

using LinearAlgebra
using DelimitedFiles
using LaTeXStrings

#################
### orbit diagram
#################

b = 0.3;
a = 1.0;
a_values = 1:0.001:1.4;
ds = Systems.henon([0.0, 0.0]; a=a, b=b)
i = 1
ics = [rand() for m in 1:10]
n = 500
Ttr = 5000

# computation
output = orbitdiagram(ds, i, 1, a_values; n = n, Ttr = Ttr)

# save data
f_name = "./tests/henon_bif_a=1.0-1.4_da=0.001_n=500_ttrans=5000.dat"
writedlm(f_name,output)
# load data
f_name = "./tests/henon_bif_a=1.0-1.4_da=0.001_n=500_ttrans=5000.dat"
output = readdlm(f_name)

# plot
pl = plot()
for (j, p) in enumerate(a_values)
    scatter!(pl, p .* ones(length(output[j,:])), output[j,:], lw = 0, ms = 0.5, color = "black")
end
plot!(pl, xlabel=L"a", ylabel=L"x", legend=nothing, xlim=[1,1.4], xticks=1:0.1:1.4, xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300)

#####################
### lyapunov exponent
#####################

sim_lyapunov_eponent = zeros(length(a_values))
for (i,a) in enumerate(a_values)
    henon = Systems.henon([0.0, 0.0]; a=a, b=b)
    @show a
    λ = lyapunov(henon, 10000; d0 = 1e-7, threshold = 1e-4, Ttr = 500)
    sim_lyapunov_eponent[i] = λ
end

# save data
f_name = "./tests/henon_lyapunov_a=1.0-1.4_da=0.001_t=10^4_ttrans=500.dat"
writedlm(f_name, sim_lyapunov_eponent)

# load data
f_name = "./tests/henon_lyapunov_a=1.0-1.4_da=0.001_t=10^4_ttrans=500.dat"
sim_lyapunov_eponent = readdlm(f_name)

pl = plot()
plot!(pl, a_values, zeros(length(a_values)), linestyle=:dash, color="black", label=nothing)
plot!(pl, a_values, sim_lyapunov_eponent, label=nothing, lw=2)
plot!(pl, xlabel=L"a", ylabel=L"\lambda", xlim=[1,1.4], xticks=1:0.1:1.4, xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300)

####################
### network measures
####################

b = 0.3;
a_values = 1:0.001:1.4;
a_values = 1.21:0.0001:1.23;
traj_length = 30000;
trans = 1000;
grid_size = 20;
ensemble = 1000;
N_steps = 10000;

sim_entropy_measures = zeros(length(a_values))
sim_lyapunov_measures = zeros(length(a_values))
theor_entropy_measures = zeros(length(a_values))
theor_lyapunov_measures = zeros(length(a_values))

#a = 1.4
#i = 1
for (i,a) in enumerate(a_values)
    system = Systems.henon([0.0, 0.0]; a=a, b=b)
    @show a
    timeseries = trajectory(system, traj_length, [0, 0]; Ttr=trans)
    discrete_timeseries, vertex_names = timeseries_to_grid(timeseries, grid_size);
    stn,retcode = create_stn(discrete_timeseries, vertex_names)
    sim_entropy_measures[i], sim_lyapunov_measures[i] = network_measures(stn, ensemble, N_steps)
    P = prob_matrix(stn);
    theor_entropy_measures[i], theor_lyapunov_measures[i] = network_measures(P)
end

data = hcat(a_values, sim_entropy_measures, sim_lyapunov_measures, theor_entropy_measures, theor_lyapunov_measures)

#save data
f_name = "./tests/henon_measures_a=1.0-1.4_da=0.001_t=10^5_ens=1000_tmax=30000_ttrans=1000.dat"
f_name = "./tests/henon_measures_a=1.21-1.23_da=0.0001_t=10^5_ens=1000_tmax=30000_ttrans=1000.dat"
writedlm(f_name,data)
# load data
f_name = "./tests/henon_measures_a=1.0-1.4_da=0.001_t=10^5_ens=1000_tmax=30000_ttrans=1000.dat"
f_name = "./tests/henon_measures_a=1.21-1.23_da=0.0001_t=10^5_ens=1000_tmax=30000_ttrans=1000.dat"
data = readdlm(f_name)

pl = scatter(data[:,1], data[:,2], label = "Random walk", lw=2, color="gray", alpha=0.5)
plot!(pl, data[:,1], data[:,4], label = "Analytic", lw=2, color="black")
plot!(pl, xlabel=L"a", ylabel=L"S", xguidefontsize=22, yguidefontsize=22, xlim=[1,1.4], xticks=1:0.1:1.4, tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300)

pl = scatter(data[:,1], data[:,3], label = "Random walk", lw=2, color="gray", alpha=0.5)
plot!(pl, data[:,1], data[:,5], label = "Analytic", lw=2, color="black")
plot!(pl, xlabel=L"a", ylabel=L"\Lambda", xguidefontsize=22, yguidefontsize=22, xlim=[1,1.4], xticks=1:0.1:1.4, tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300)

#############
### full plot
#############

a_values = 1:0.001:1.4;
a_values_zoom = 1.21:0.0001:1.23;
# load bif data
f_name = "./tests/henon_bif_a=1.0-1.4_da=0.001_n=500_ttrans=5000.dat"
output = readdlm(f_name)
# load measures data
f_name = "./tests/henon_measures_a=1.0-1.4_da=0.001_t=10^5_ens=1000_tmax=30000_ttrans=1000.dat"
data = readdlm(f_name)
f_name = "./tests/henon_measures_a=1.21-1.23_da=0.0001_t=10^5_ens=1000_tmax=30000_ttrans=1000.dat"
data_zoom = readdlm(f_name)
# load lyaunov exponent data
f_name = "./tests/henon_lyapunov_a=1.0-1.4_da=0.001_t=10^4_ttrans=500.dat"
sim_lyapunov_eponent = readdlm(f_name)

l = @layout [a{0.4h}; b{0.2h}; c{0.2h}; b{0.2h}]

p1 = plot()
for (j, p) in enumerate(a_values)
    scatter!(p1, p .* ones(length(output[j,:])), output[j,:], lw = 0, ms = 0.5, color = "black")
end
scatter!(p1, [1.2], [1.3], label=nothing, lw=0, ms=10, markerstrokewidth=0, color="orange")
scatter!(p1, [1.2265], [1.3], label=nothing, lw=0, ms=10, markerstrokewidth=0, color="red")
scatter!(p1, [1.24], [1.3], label=nothing, lw=0, ms=10, markerstrokewidth=0, color="blue")
scatter!(p1, [1.27], [1.3], label=nothing, lw=0, ms=10, markerstrokewidth=0, color="green")
plot!(p1, ylabel=L"x_n", legend=nothing, xlim=[1,1.4], xticks=1:0.1:1.4, xguidefontsize=22, yguidefontsize=22, tickfontsize=14, xtickfontsize=1, lw=2, legendfontsize=16,

dpi=300)

p2 = plot()
plot!(p2, a_values, zeros(length(a_values)), linestyle=:dash, color="black", label=nothing)
plot!(p2, a_values, sim_lyapunov_eponent, label=nothing, color="black", lw=2)
plot!(p2, ylabel=L"\lambda", xlim=[1,1.4], xticks=1:0.1:1.4, yticks=-0.5:0.5:0.5, xguidefontsize=22, yguidefontsize=22, tickfontsize=14, xtickfontsize=1, lw=2, fontfamily="serif", legendfontsize=16,
dpi=300)

p3 = plot()
scatter!(p3, data[:,1], data[:,2], label = "Random walk", lw=0, ms=6, markerstrokewidth=0, color="gray", alpha=0.5)
scatter!(p3, data_zoom[:,1], data_zoom[:,2], label=nothing, lw=0, ms=6, markerstrokewidth=0, color="red", alpha=0.2)
plot!(p3, data[:,1], data[:,4], label = "Analytic", lw=2, color="black")
plot!(p3, data_zoom[:,1], data_zoom[:,4], label=nothing, lw=2, color="red")
plot!(p3, ylabel=L"S", xguidefontsize=22, yguidefontsize=22, xlim=[1,1.4], xticks=1:0.1:1.4, yticks=0:0.5:1., ylim=[0,1],tickfontsize=14, xtickfontsize=1, lw=2, fontfamily="serif", legendfontsize=16, legend_position=:topleft,
dpi=300)

p4 = plot()
scatter!(p4, data[:,1], data[:,3], label = "Random walk", lw=0, ms=6, markerstrokewidth=0, color="gray", alpha=0.5)
scatter!(p4, data_zoom[:,1], data_zoom[:,3], label=nothing, lw=0, ms=6, markerstrokewidth=0, color="red", alpha=0.2)
plot!(p4, data[:,1], data[:,5], label = "Analytic", lw=2, color="black")
plot!(p4, data_zoom[:,1], data_zoom[:,5], label=nothing, lw=2, color="red")
plot!(p4, xlabel=L"a", ylabel=L"\Lambda", xguidefontsize=22, yguidefontsize=22, xlim=[1,1.4], xticks=1:0.1:1.4, yticks=0:0.5:1, tickfontsize=14, lw=2, fontfamily="serif", legendfontsize=16,
dpi=300)

p = plot(p1, p2, p3, p4, size=(1000,1000), layout = l);
savefig(p, "./tests/henon_measures_a=1.0-1.4_da=0.001.png")

# correlation plots
pl = plot()
scatter!(pl, sim_lyapunov_eponent, data[:,4], lw=0, ms=6, markerstrokewidth=0, color="purple", alpha=0.5)
plot!(pl, legend=nothing, xlabel=L"\lambda", ylabel=L"S", xlim=[0,0.4], xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="serif", legendfontsize=16, dpi=300)
savefig(pl, "./tests/henon_corr_lyapexp-S.png")

pl = plot()
scatter!(pl, sim_lyapunov_eponent, data[:,5], lw=0, ms=6, markerstrokewidth=0, color="cyan", alpha=0.5)
plot!(pl, legend=nothing, xlabel=L"\lambda", ylabel=L"\Lambda", xlim=[0.,0.4],xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="serif", legendfontsize=16, dpi=300)
savefig(pl, "./tests/henon_corr_lyapexp-lyapmes.png")


# scaling plot
a_c = 1.2266
pl = plot()
scatter!(pl, abs.(data_zoom[:,1] .- a_c)/a_c, data_zoom[:,5] ./ data_zoom[:,4], label=L"a_c=%$(a_c)", lw=0, ms=6, markerstrokewidth=0, color="red", alpha=0.5)
plot!(pl, xaxis=:log, yaxis=:log, xlabel=L"(a-a_c)/a_c", ylabel=L"\Lambda/S", xlim=[1e-5,1e-1], xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="sans-serif", legendfontsize=16, dpi=300)
savefig(pl, "./tests/henon_scaling_ac=1.2266.png")