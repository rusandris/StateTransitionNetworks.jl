using DynamicalSystemsBase
using StateTransitionNetworks
using Plots, LaTeXStrings
using DelimitedFiles

#------------------------define duffing system with parameter drift------------------------

init_drift_pars = [1.0, 0.2, 0.4, 0.0005]
init_nodrift_pars = [1.0, 0.2, 0.4, 0.0]

function duffing_drift(u0 = [0.1, 0.25]; init_sys_pars=init_drift_pars)
    ω = init_sys_pars[1]
    β = init_sys_pars[2]
    ε0 = init_sys_pars[3]
    α = init_sys_pars[4]
    return CoupledODEs(duffing_drift_rule, u0, [ω, β, ε0, α])
end

@inbounds function duffing_drift_rule(x, p, t)
    ω, β, ε0, α = p
    dx1 = x[2]
    dx2 = (ε0+α*t)*cos(ω*t) + x[1] - x[1]^3 - 2β * x[2]
    return SVector(dx1, dx2)
end

# create continuous time system
ds = duffing_drift(; init_sys_pars=init_nodrift_pars)

# create map from continuous system
#duffing_map = StroboscopicMap(ds,2π)

# parameters for testing
eps_range = [0.4:0.01:0.8;] #[0.4:0.001:0.8;]
N_ens = 80^2 # 200^2 #80^2
init_ensemble_grid_edges=[(-1.5,1.5),(-1.5,1.5)]
Ttr = 400

#---------------------run ensemble symbolics, no drift--------------------------
init_sys_pars = init_nodrift_pars
results_dir = "./results/"; mkpath(results_dir)

#----------grid partition------------
symbolization_pars_grid = (symb_type=:grid,grid_size=20,grid_edges=[-1.5,-1.5,1.5,1.5])

#order of symbols 
order0 = 1
Ss_ens_order0, Λs_ens_order0, _, _ = ensemble_measures(eps_range,N_ens,duffing_drift,symbolization_pars_grid,order0;
               init_ensemble_grid_edges=init_ensemble_grid_edges,Ttr=Ttr,
               init_sys_pars=init_sys_pars,p_index=3,Δt=2π/init_drift_pars[1])
               
writedlm(results_dir*"duffing_nodrift_measuresSL-n_omega=1_eps0=0.4_beta=0.2_alpha=$(init_sys_pars[end])_Nens=$(N_ens)_Ttr=$(Ttr)_grid=$(symbolization_pars_grid.grid_size)_order=$(order0)_ergodic.txt",hcat(Ss_ens_order0, Λs_ens_order0))

#order of symbols 
order1 = 2
Ss_ens_order1, Λs_ens_order1, _, _ = ensemble_measures(eps_range,N_ens,duffing_drift,symbolization_pars_grid,order1;
               init_ensemble_grid_edges=init_ensemble_grid_edges,Ttr=Ttr,
               init_sys_pars=init_sys_pars,p_index=3,Δt=2π/init_drift_pars[1])

writedlm(results_dir*"duffing_nodrift_measuresSL-n_omega=1_eps0=0.4_beta=0.2_alpha=$(init_sys_pars[end])_Nens=$(N_ens)_Ttr=$(Ttr)_grid=$(symbolization_pars_grid.grid_size)_order=$(order1)_ergodic.txt",hcat(Ss_ens_order1, Λs_ens_order1))

#order of symbols 
order2 = 3            
Ss_ens_order2, Λs_ens_order2, _, _ = ensemble_measures(eps_range,N_ens,duffing_drift,symbolization_pars_grid,order2;
               init_ensemble_grid_edges=init_ensemble_grid_edges,Ttr=Ttr,
               init_sys_pars=init_sys_pars,p_index=3,Δt=2π/init_drift_pars[1])
               
writedlm(results_dir*"duffing_nodrift_measuresSL-n_omega=1_eps0=0.4_beta=0.2_alpha=$(init_sys_pars[end])_Nens=$(N_ens)_Ttr=$(Ttr)_grid=$(symbolization_pars_grid.grid_size)_order=$(order2)_ergodic.txt",hcat(Ss_ens_order2, Λs_ens_order2))               

#---------------------plot comparing entropies and lyapunov exponents--------------------------------
pl_grid_S = plot();

data_lyaps = readdlm("./duffing_drift_lyapunov-eps0_omega=1_beta=0.2_alpha=0_Nens=10^3_Ttr=200.dat")
plot!(pl_grid_S,data_lyaps[1:end,1],data_lyaps[1:end,2],lc=:gray10,lw=1,alpha=0.8,label=L"\lambda_\mathrm{ens}",dpi=300);

#plot!(pl_grid_S,eps_range,Ss_ens_order1;lc=:lightgoldenrod1,lw=2,alpha=0.8,label=L"S^{grid}_{ens}(%$(order0))");
plot!(pl_grid_S,eps_range,Ss_ens_order1;lc=:orange,lw=2,alpha=0.8,label=L"S^{grid}_{ens}(%$(order1))");
#plot!(pl_grid_S,eps_range,Ss_ens_order2;lc=:brown,lw=2,alpha=0.8,label=L"S^{grid}_{ens}(%$(order2))");
plot!(pl_grid_S,xlabel=L"\varepsilon_0",ylabel=L"S,\lambda", ylims=[-0.05,1.6], xguidefontsize=22, 
	yguidefontsize=22, tickfontsize=14, legendfontsize=16);

pl_grid_L = plot();
#plot!(pl_grid_L,eps_range,Λs_ens_order0;lc=:salmon,label=L"\Lambda^{grid}_{ens}(%$(order0))");
plot!(pl_grid_L,eps_range,Λs_ens_order1;lc=:red,label=L"\Lambda^{grid}_{ens}(%$(order1))");
#plot!(pl_grid_L,eps_range,Λs_ens_order2;lc=:purple,label=L"\Lambda^{grid}_{ens}(%$(order2))");
plot!(pl_grid_L,xlabel=L"\varepsilon_0",ylabel=L"\Lambda", ylims=[-0.05,1.], xguidefontsize=22, 
	yguidefontsize=22, tickfontsize=14, legendfontsize=16);

pm1 = plot(pl_grid_S, pl_grid_L, size=(800,600), dpi=600, layout=grid(2, 1, heights=[0.5, 0.5]), legend=:bottomright)

figs_dir = "./figs/";mkpath(figs_dir)
savefig(pm1,figs_dir*"duffing_nodrift_S-instlyap-n_omega=1_eps0=0.4_beta=0.2_alpha=$(init_sys_pars[end])_Nens=$(N_ens)_Ttr=$(Ttr)_grid=$(symbolization_pars_grid.grid_size)_ergodic.pdf")


#---------------------run ensemble symbolics, with drift--------------------------
#order of symbols 
order = 2
n_steps = 130

#----------grid partition------------
symbolization_pars_grid = (symb_type=:grid,grid_size = 2^3,grid_edges=[-1.5,-1.5,1.5,1.5])

Ss_ens_drift, Λs_ens_drift, _, _ = ensemble_measures(n_steps,N_ens,duffing_drift,symbolization_pars_grid,order;
               init_ensemble_grid_edges=init_ensemble_grid_edges,Ttr=Ttr,
               auto_sys_pars=init_nodrift_pars,nonauto_sys_pars=init_drift_pars,p_index=3,Δt=2π/init_drift_pars[1])

writedlm(results_dir*"duffing_drift_measuresSL-n_omega=1_eps0=0.4_beta=0.2_alpha=$(init_sys_pars[end])_Nens=$(N_ens)_Ttr=$(Ttr)_grid=$(symbolization_pars_grid.grid_size)_order=$(order2)_ergodic.txt",hcat(Ss_ens_drift, Λs_ens_drift))  

#------------plot drift results------------
steps = 0:n_steps
drift_param_vals = init_drift_pars[3] .+ init_drift_pars[4] .* steps .* 2π


pl_grid_S = plot();
data_lyaps = readdlm("./duffing_drift_instlyap-n_omega=1_eps0=0.4_beta=0.2_alpha=0.0005_Nens=40000_Ttr=200_maxrescaletime=2.dat")
plot!(pl_grid_S,data_lyaps[1:end,1],data_lyaps[1:end,2],lc=:gray10,lw=1,alpha=0.8,label=L"\lambda_\mathrm{ens}",dpi=300);

plot!(pl_grid_S, [0,130], [0,0], label=nothing, ls=:dash, color=:black, lw=1., alpha=0.8);
plot!(pl_grid_S, xlabel=L"n",ylabel=L"S, \lambda", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, legendfontsize=16)
plot!(pl_grid_S, ylims=[-0.1, 1.6], xguidefontsize=22, yguidefontsize=22, tickfontsize=14, legendfontsize=16)
plot!(twiny(), drift_param_vals, Ss_ens_drift, framestyle=:semi, label=false, xlabel=L"\varepsilon", yticks=false, color=:red, lw=2, alpha=0.);

plot!(pl_grid_S,steps,Ss_ens_drift;lc=:orange,lw=2,alpha=0.8,label=L"S^{grid}_{ens}(%$(order))");

pl_grid_L = plot();
data_lyaps = readdlm("./duffing_drift_instlyap-n_omega=1_eps0=0.4_beta=0.2_alpha=0.0005_Nens=40000_Ttr=200_maxrescaletime=2.dat")
#plot!(pl_grid_op,data_lyaps[1:end,1],data_lyaps[1:end,2],lc=:gray10,lw=1,alpha=0.8,label=L"\lambda_\mathrm{ens}",dpi=300);

plot!(pl_grid_L, [0,130], [0,0], label=nothing, ls=:dash, color=:black, lw=1., alpha=0.8);
plot!(pl_grid_L, xlabel=L"n",ylabel=L"S, \lambda", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, legendfontsize=16)
plot!(pl_grid_L, ylims=[-0.1, 1.6], xguidefontsize=22, yguidefontsize=22, tickfontsize=14, legendfontsize=16)
plot!(twiny(), drift_param_vals, Ss_ens_drift, framestyle=:semi, label=false, xlabel=L"\varepsilon", yticks=false, color=:red, lw=2, alpha=0.);

plot!(pl_grid_L,steps,Λs_ens_drift;lc=:red,lw=2,alpha=0.8,label=L"\Lambda^{grid}_{ens}(%$(order))");

pm2 = plot(pl_grid_S, pl_grid_L, size=(800,600), dpi=600, layout=grid(2, 1, heights=[0.5, 0.5]), legend=:bottomright)
            
savefig(pm2,figs_dir*"duffing_drift_S-instlyap-n_omega=1_eps0=0.4_beta=0.2_alpha=$(init_drift_pars[end])_Nens=$(N_ens)_Ttr=$(Ttr)__grid=$(symbolization_pars_grid.grid_size)_order=$(order)_ergodic.pdf")
