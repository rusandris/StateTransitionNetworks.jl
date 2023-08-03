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
plane = (1,15.0);
u0 = rand(3)
rho_values = 180.:0.001:182.;
T = 20000;
ds = PredefinedDynamicalSystems.lorenz(u0)
Ttr = 19990

# computation
output = produce_orbitdiagram(ds, plane, 3, 2, rho_values; tfinal=T, Ttr=Ttr, direction=1, printparams=true, rootkw=(xrtol=1e-10, atol=1e-10))

# save data
f_name = "./tests/lorenz_bif_rho=180-182_drho=0.001_tmax=20000_ttrans=19900.dat"
#writedlm(f_name,output)
serialize(f_name, output)
# load data
f_name = "./tests/lorenz_bif_rho=180-182_drho=0.001_tmax=20000_ttrans=19900.dat"
#output = readdlm(f_name)
output = deserialize(f_name)

# plot
pl = plot()
for (j, p) in enumerate(rho_values)
    mask = (output[j,:] .!= "")
    scatter!(pl, p .* ones(Int(sum(mask))), output[j,mask] .- 180, lw = 0, ms = 0.5, color = "black")
end
plot!(pl, xlabel=L"\rho", ylabel=L"z_n-180", legend=nothing, xlim=[180,182], ylim=[3.5,5.5], xticks=180:0.5:182, xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="serif", legendfontsize=16,
dpi=300);
savefig(pl, "./tests/lorenz_bif.png")


#####################
### lyapunov exponent
#####################

sim_lyapunov_eponent = zeros(length(rho_values))
for (i,ρ) in enumerate(rho_values)
    ds = PredefinedDynamicalSystems.lorenz(; ρ=ρ)
    @show ρ
    λ = lyapunov(ds, 10000; d0=1e-7, Ttr=5000)
    sim_lyapunov_eponent[i] = λ
end

# save data
f_name = "./tests/lorenz_lyapunov_rho=180-182_drho=0.001_t=10^4_ttrans=5000.dat"
writedlm(f_name, sim_lyapunov_eponent)

# load data
f_name = "./tests/lorenz_lyapunov_rho=180-182_drho=0.001_t=10^4_ttrans=5000.dat"
sim_lyapunov_eponent = readdlm(f_name)

pl = plot()
plot!(pl, rho_values, zeros(length(rho_values)), linestyle=:dash, color="black", label=nothing)
plot!(pl, rho_values, sim_lyapunov_eponent, label=nothing, lw=2)
plot!(pl, xlabel=L"\rho", ylabel=L"\lambda", xlim=[180,182], xticks=180:0.5:182, xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="serif", legendfontsize=16,
dpi=300)

####################
### network measures
####################
Δt = 0.001;
plane = (1,15.0);
grid = 20;
u0 = rand(3)
rho_values = 180.:0.001:182.;
rho_values = 181.66:0.0001:181.72;
T = 100000;
trans = 5000
ensemble = 1000;
N_steps = 10000;
grid_size = 20;

sim_entropy_measures = zeros(length(rho_values))
sim_lyapunov_measures = zeros(length(rho_values))
theor_entropy_measures = zeros(length(rho_values))
theor_lyapunov_measures = zeros(length(rho_values))

ρ = 180.1;
i = 1
stn=nothing
for (i,ρ) in enumerate(rho_values)
        @show ρ
        while true # STN must be 'healthy'
            u0 = rand(Float64,3).*50 .-25;
            system = PredefinedDynamicalSystems.lorenz(u0; ρ=ρ);
            timeseries,  = trajectory(system, T; Δt=Δt, Ttr=trans);
            psection = DynamicalSystemsBase.poincaresos(timeseries, plane; direction=+1, save_idxs=[2,3]);
            d_traj, v_names = timeseries_to_grid(psection, grid);
            stn, ret_code_stn = create_stn(d_traj, v_names);
            if ret_code_stn == :Success
                @show stn
                break
            end
        end
    sim_entropy_measures[i], sim_lyapunov_measures[i] = network_measures(stn, ensemble, N_steps)
    P = prob_matrix(stn);
    theor_entropy_measures[i], theor_lyapunov_measures[i] = network_measures(P)
end

data = hcat(rho_values, sim_entropy_measures, sim_lyapunov_measures, theor_entropy_measures, theor_lyapunov_measures)
f_name = "./tests/lorenz_measures_rho=180-182_drho=0.001_ens=1000_tmax=10^5_ttrans=5000.dat"
f_name = "./tests/lorenz_measures_rho=181.65-181.73_drho=0.001_ens=1000_tmax=10^5_ttrans=5000.dat"
f_name = "./tests/lorenz_measures_rho=181.66-181.72_drho=0.0001_ens=1000_tmax=10^5_ttrans=5000.dat"
writedlm(f_name,data)

f_name = "./tests/lorenz_measures_rho=180-182_drho=0.001_ens=1000_tmax=10^5_ttrans=5000.dat"
f_name = "./tests/lorenz_measures_rho=181.65-181.73_drho=0.001_ens=1000_tmax=10^5_ttrans=5000.dat"
f_name = "./tests/lorenz_measures_rho=181.66-181.72_drho=0.0001_ens=1000_tmax=10^5_ttrans=5000.dat"
data = readdlm(f_name)

pl = scatter(data[:,1], data[:,2], label = "Random walk", lw=2, color="gray", alpha=0.5)
plot!(pl, data[:,1], data[:,4], label = "Analytic", lw=1, color="black")
plot!(pl, xlabel=L"\rho", ylabel=L"S", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300)

pl = scatter(data[:,1], data[:,3], label = "Random walk", lw=1, ms=5, color="gray", alpha=0.5)
plot!(pl, data[:,1], data[:,5], label = "Analytic", lw=1, color="black", alpha=0.5)
plot!(pl, xlabel=L"\rho", ylabel=L"\Lambda", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300)

#############
### full plot
#############

rho_values = 180:0.001:182;
rho_values_zoom = 181.65:0.001:181.73;
# load bif data
f_name = "./tests/lorenz_bif_rho=180-182_drho=0.001_tmax=20000_ttrans=19900.dat"
output = deserialize(f_name)
# load measures data
f_name = "./tests/lorenz_measures_rho=180-182_drho=0.001_ens=1000_tmax=10^5_ttrans=5000.dat"
data = readdlm(f_name)
f_name = "./tests/lorenz_measures_rho=181.65-181.73_drho=0.001_ens=1000_tmax=10^5_ttrans=5000.dat"
f_name = "./tests/lorenz_measures_rho=181.66-181.72_drho=0.0001_ens=1000_tmax=10^5_ttrans=5000.dat"
data_zoom = readdlm(f_name)
# load lyaunov exponent data
f_name = "./tests/lorenz_lyapunov_rho=180-182_drho=0.001_t=10^4_ttrans=5000.dat"
sim_lyapunov_eponent = readdlm(f_name)

l = @layout [a{0.4h}; b{0.2h}; c{0.2h}; b{0.2h}]

p1 = plot();
for (j, p) in enumerate(rho_values)
    mask = (output[j] .!= "")
    scatter!(p1, p .* ones(Int(sum(mask))), output[j][mask] .- 180., lw=0, ms=0.2, color="black", alpha=0.1)
end
scatter!(p1, [180.10], [5.3], label=nothing, lw=0, ms=10, markerstrokewidth=0, color="orange");
scatter!(p1, [180.70], [5.3], label=nothing, lw=0, ms=10, markerstrokewidth=0, color="red");
scatter!(p1, [181.1], [5.3], label=nothing, lw=0, ms=10, markerstrokewidth=0, color="blue");
scatter!(p1, [180.78], [5.3], label=nothing, lw=0, ms=10, markerstrokewidth=0, color="green");
plot!(p1, ylabel=L"z_n-180", legend=nothing, xlim=[180,182], xticks=180:0.5:182, ylim=[3.5,5.5], xguidefontsize=22, yguidefontsize=22, tickfontsize=14, xtickfontsize=1, lw=2, legendfontsize=16, dpi=300);

rho_values = 180:0.001:182;
p2 = plot();
plot!(p2, rho_values, zeros(length(rho_values)), linestyle=:dash, color="black", label=nothing);
plot!(p2, rho_values, sim_lyapunov_eponent, label=nothing, color="black", lw=2);
plot!(p2, ylabel=L"\lambda", xlim=[180,182], xticks=180:0.5:182, yticks=0:1.:2., xguidefontsize=22, yguidefontsize=22, tickfontsize=14, xtickfontsize=1, lw=2, fontfamily="serif", legendfontsize=16, dpi=300);

p3 = plot();
scatter!(p3, data[:,1], data[:,2], label = "Random walk", lw=0, ms=6, markerstrokewidth=0, color="gray", alpha=0.5);
scatter!(p3, data_zoom[:,1], data_zoom[:,2], label=nothing, lw=0, ms=6, markerstrokewidth=0, color="red", alpha=0.2);
plot!(p3, data[:,1], data[:,4], label = "Analytic", lw=2, color="black");
plot!(p3, data_zoom[:,1], data_zoom[:,4], label=nothing, lw=2, color="red");
plot!(p3, ylabel=L"S", xguidefontsize=22, yguidefontsize=22, xlim=[180,182], xticks=180:0.5:182, yticks=0:1.:2., ylim=[0,2],tickfontsize=14, xtickfontsize=1, lw=2, fontfamily="serif", legendfontsize=16, legend_position=:bottomleft, dpi=300);

p4 = plot();
scatter!(p4, data[:,1], data[:,3], label = "Random walk", lw=0, ms=6, markerstrokewidth=0, color="gray", alpha=0.5);
scatter!(p4, data_zoom[:,1], data_zoom[:,3], label=nothing, lw=0, ms=6, markerstrokewidth=0, color="red", alpha=0.2);
plot!(p4, data[:,1], data[:,5], label = "Analytic", lw=2, color="black");
plot!(p4, data_zoom[:,1], data_zoom[:,5], label=nothing, lw=2, color="red");
plot!(p4, xlabel=L"\rho", ylabel=L"\Lambda", xguidefontsize=22, yguidefontsize=22, xlim=[180,182], xticks=180:0.5:182, yticks=0:1.:4, tickfontsize=14, lw=2, fontfamily="serif", legendfontsize=16, dpi=300);

p = plot(p1, p2, p3, p4, size=(1000,1000), layout = l);
savefig(p, "./tests/lorenz_measures_rho=180-182_drho=0.001.png")

# correlation plots
pl = plot()
scatter!(pl, sim_lyapunov_eponent, data[:,4], lw=0, ms=6, markerstrokewidth=0, color="purple", alpha=0.5)
plot!(pl, legend=nothing, xlabel=L"\lambda", ylabel=L"S", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="sans-serif", legendfontsize=16, dpi=300)
savefig(pl, "./tests/lorenz_corr_lyapexp-S.png")

pl = plot()
scatter!(pl, sim_lyapunov_eponent, data[:,5], lw=0, ms=6, markerstrokewidth=0, color="cyan", alpha=0.5)
plot!(pl, legend=nothing, xlabel=L"\lambda", ylabel=L"\Lambda", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="serif", legendfontsize=16, dpi=300)
savefig(pl, "./tests/lorenz_corr_lyapexp-lyapmes.png")


# scaling plot
ρ_c = 181.669
pl = plot()
scatter!(pl, (data_zoom[:,1] .- ρ_c)/ρ_c, data_zoom[:,5] ./ data_zoom[:,4], label=L"\rho_c=%$(ρ_c)", lw=0, ms=6, markerstrokewidth=0, color="red", alpha=0.5)
plot!(pl, xaxis=:log, yaxis=:log, xlabel=L"(\rho-\rho_c)/\rho_c", ylabel=L"\Lambda/S", xlim=[1e-6,1e-3], xguidefontsize=22, yguidefontsize=22, tickfontsize=14, lw=2, fontfamily="sans-serif", legendfontsize=16, dpi=300)
savefig(pl, "./tests/lorenz_scaling_rhoc=181.669.png")