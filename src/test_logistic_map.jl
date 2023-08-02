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

function transition_matrix(f, grid, grid_edges)
    P = zeros(grid,grid)
    x_min,x_max = grid_edges 
    dx = 0.5*(x_max-x_min)/grid
    x = collect(range(grid_edges[1],grid_edges[2], grid))
    for i in 1:grid
        P[i,Int(f(i))]
    end
end

function bit_number_measures(P::AbstractMatrix)
	λ, X = eigen(transpose(Matrix(P)))
	
	if real(λ[end]) ≈ 1
	   x = transpose(X[:,end]./sum(X[:,end]))
	else
	   return -1, :StochasticMatrixError
	end

	l = -log.(x)
	replace!(l, Inf=>0.0)
	entropy = sum((x .* l))
    #variance = sum(x .* ((l .- entropy).^2))
    variance =  sum(x .* l .* l) - entropy^2
	return real(entropy), real(variance), :Success
end

###############################
### Measures for a single value
###############################

grid_size = 20;
r = 4.;
ds = Systems.logistic(0.4; r = r)
λ = lyapunov(ds, 10000; Ttr = 500)
timeseries = trajectory(ds, 10^8, 0.4; Ttr=1000);
ys = zeros(length(timeseries))
timeseries = hcat(timeseries, ys)
discrete_timeseries, vertex_names = timeseries_to_grid(timeseries, grid_size; grid_edges = [0., 0., 1., 1.]);
@time stn,retcode = create_stn(discrete_timeseries, vertex_names);
P = prob_matrix(stn);
@time S, Λ = network_measures(P)
@time C1, C2 = bit_number_measures(P) 

####################
### network measures
####################

r_values = 3.995:0.001:4.;
traj_length = 10^8;
trans = 1000;
grid_size = 20;
ensemble = 1000;
N_steps = 10000;

sim_entropy_measures = zeros(length(r_values))
sim_lyapunov_measures = zeros(length(r_values))
theor_entropy_measures = zeros(length(r_values))
theor_lyapunov_measures = zeros(length(r_values))
C1_measures = zeros(length(r_values))
C2_measures = zeros(length(r_values))

for (i,r) in enumerate(r_values)
    system = Systems.logistic(0.4; r=r)
    @show r
    timeseries = trajectory(system, traj_length, 0.4; Ttr=trans)
    ys = zeros(length(timeseries))
    timeseries = hcat(timeseries, ys)
    discrete_timeseries, vertex_names = timeseries_to_grid(timeseries, grid_size; grid_edges = [0., 0., 1., 1.]);
    stn,retcode = create_stn(discrete_timeseries, vertex_names)
    #sim_entropy_measures[i], sim_lyapunov_measures[i] = network_measures(stn, ensemble, N_steps)
    P = prob_matrix(stn);
    theor_entropy_measures[i], theor_lyapunov_measures[i] = network_measures(P)
    #C1_measures[i], C2_measures[i] = bit_number_measures(P)
end

data = hcat(r_values, sim_entropy_measures, sim_lyapunov_measures, theor_entropy_measures, theor_lyapunov_measures, C1_measures, C2_measures)

#save data
f_name = "./tests/logistic_measures_r=3.5-4_da=0.001_ens=1000_tmax=100000_ttrans=1000_grid=200.dat"
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20.dat"
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=10^2_ttrans=1000_grid=200.dat"
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=10^3_ttrans=1000_grid=200.dat"
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=10^4_ttrans=1000_grid=200.dat"
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=10^5_ttrans=1000_grid=200.dat"
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=200.dat"
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=10^7_ttrans=1000_grid=200.dat"
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=1000.dat"
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=1000000_ttrans=1000_grid=2000.dat"
writedlm(f_name,data)
# load data
f_name = "./tests/logistic_measures_r=0-4_da=0.01_ens=1000_tmax=100000_ttrans=1000_grid=200.dat"
data = readdlm(f_name)

# plot bit number measures C2
pl = plot()
plot!(pl, data[:,1], data[:,5], label = L"\Lambda", lw=2, alpha=0.5, color="red")
plot!(pl, data[:,1], data[:,7], label = L"C_2", lw=2, alpha=0.5, color="orange")
plot!(pl, xlabel=L"r", ylabel=L"\Lambda, C_2", xguidefontsize=22, yguidefontsize=22, xlim=[3.5,4], tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300)
plot!(ylim=[0,1], xlim=[3.6,3.7])
savefig(pl, "./tests/logistic_var_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20.svg")
savefig(pl, "./tests/logistic_var_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20.pdf")

# plot bit number measures C1
pl = plot()
plot!(pl, data[:,1], data[:,4], label = L"S", lw=2, alpha=0.5, color="black")
plot!(pl, data[:,1], data[:,6], label = L"C_1", lw=2, alpha=0.5, color="brown")
plot!(pl, xlabel=L"r", ylabel=L"S, C_1", xguidefontsize=22, yguidefontsize=22, xlim=[3.5,4], tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300)
savefig(pl, "./tests/logistic_entr_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20.svg")
savefig(pl, "./tests/logistic_entr_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20.pdf")

pl = plot()
#scatter!(pl, data[:,1], data[:,2], label = "Random walk", lw=2, color="gray", alpha=0.5)
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=1000.dat"
data = readdlm(f_name)
plot!(pl, data[:,1], data[:,4], label = "grid=1000", lw=2, alpha=0.5, color="black")
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=200.dat"
data = readdlm(f_name)
plot!(pl, data[:,1], data[:,4], label = "grid=200", lw=2, alpha=0.5, color="red")
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20.dat"
data = readdlm(f_name)
plot!(pl, data[:,1], data[:,4], label = "grid=20", lw=2, alpha=0.5, color="blue")
plot!(pl, xlabel=L"r", ylabel=L"S", xguidefontsize=22, yguidefontsize=22, xlim=[3.5,4], tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300)
savefig(pl, "./tests/logistic_entr_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20-200-1000.svg")
savefig(pl, "./tests/logistic_entr_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20-200-1000.pdf")


# lyapunov measure
pl = plot()
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=1000.dat"
data = readdlm(f_name)
plot!(pl, data[:,1], data[:,5].+0.01, label = "grid=1000", lw=2, alpha=0.5, color="black")
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=200.dat"
data = readdlm(f_name)
plot!(pl, data[:,1], data[:,5].+0.01, label = "grid=200", lw=2, alpha=0.5, color="red")
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20.dat"
data = readdlm(f_name)
plot!(pl, data[:,1], data[:,5].+0.01, label = "grid=20", lw=2, alpha=0.5, color="blue")
plot!(pl, xlabel=L"r", ylabel=L"\Lambda", xguidefontsize=22, yguidefontsize=22, xlim=[3.5,4], tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300, yscale=:linear)
plot!(ylim=[0,2])
savefig(pl, "./tests/logistic_lyap_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20-200-1000.svg")
savefig(pl, "./tests/logistic_lyap_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20-200-1000.pdf")

# compare gris in r\in [3.6,3.7]
pl = plot()
#scatter!(pl, data[:,1], data[:,3], label = "Random walk", lw=2, color="gray", alpha=0.5)
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=1000.dat"
data = readdlm(f_name)
plot!(pl, data[:,1], data[:,5].+0.01, label = "grid=1000", lw=2, alpha=0.5, color="black")
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=200.dat"
data = readdlm(f_name)
plot!(pl, data[:,1], data[:,5].+0.01, label = "grid=200", lw=2, alpha=0.5, color="red")
f_name = "./tests/logistic_measures_r=3.5-4_da=0.0001_tmax=100000_ttrans=1000_grid=20.dat"
data = readdlm(f_name)
plot!(pl, data[:,1], data[:,5].+0.01, label = "grid=20", lw=2, alpha=0.5, color="blue")
plot!(pl, xlabel=L"r", ylabel=L"\Lambda+0.01", xguidefontsize=22, yguidefontsize=22, xlim=[3.6,3.7], tickfontsize=14, lw=2, fontfamily="Serif", legendfontsize=16,
dpi=300, yscale=:log)
savefig(pl, "./tests/logistic_lyap_r=3.6-3.7_da=0.0001_tmax=100000_ttrans=1000_grid=20-200-1000.svg")
savefig(pl, "./tests/logistic_lyap_r=3.6-3.7_da=0.0001_tmax=100000_ttrans=1000_grid=20-200-1000.pdf")