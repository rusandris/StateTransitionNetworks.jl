using StateTransitionNetworks
using Plots
using DynamicalSystems

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


guidefontsize=20
tickfontsize=15
legendfontsize=15


Δt = 0.001;
plane = (2,0.0);
grid = 20;
T = 5000;
ensemble = 10000
N_steps = 10000
L_values = []

b_vals = [0.42,0.368,0.28]
logocolors = Colors.JULIA_LOGO_COLORS
xticks_forall = [[8000,8200,8400],[2200,2600,3000],[6150,6300,6450]]

#rho_arr = [180.1]
color_arr = [:orange, :red, :green]
individual_plots = []
pl = plot()

lyaps = []
entropies = []

for (i,b) in enumerate(b_vals)
    system = PredefinedDynamicalSystems.roessler(b=b);
    timeseries,  = trajectory(system, T; Δt=Δt, Ttr=500);
    psection = DynamicalSystemsBase.poincaresos(timeseries, plane; direction=+1, save_idxs=[1,3]);
    d_traj, v_names = timeseries_to_grid(psection, grid);
    stn, ret_code = create_stn(d_traj, v_names);
    
    #=
    P = prob_matrix(P)
    entropy,lyapunov = network_measures(P)
    =#
    
    ret_code
    L, S, Λ = walk_length_distribution(stn, ensemble, N_steps)
    
    

    S_hist = mean(L)/N_steps
    Λ_hist = var(L, corrected=false)/N_steps
    
    push!(lyaps,Λ_hist)
	push!(entropies,S_hist)
	
    μ = S_hist*N_steps
    σ = sqrt(Λ_hist*N_steps)
    L_range = collect(range(μ-4σ,μ+4σ,1000))
    

	#------------------three plots together----------------

    
    pl_b = plot()
    
    ylabel = i == 1 ? L"$p(L(t))$" : ""
    yticks = i == 1 ? [0.0,5e-3,1e-2] : false
    
    histogram!(pl, L, normalize=:pdf, bins=100,
        color=color_arr[i],
        alpha=0.3,
        lw = 0.00001,
        la=0.001,
        xlabel=L"L(t)",
        ylabel=L"$p(L(t))$",
       	guidefontsize=guidefontsize,
		tickfontsize=tickfontsize,
		legendfontsize=legendfontsize,
        #label=L"\rho_{%$(i)}=%$(ρ)",
        label=L"b="*"$b",
        legend=:topleft,
        framestyle=:box,
        dpi=300,
        ylims=[0,1.2e-2],
        yticks = [0.0,5e-3,1e-2],
        xticks = [2000,5500,9000],
        margin=8Plots.mm)
        
        plot!(pl, L_range, gauss.(L_range,μ,σ),
        color=color_arr[i],
        linewidth=2,
        label="",
        framestyle=:box,
        dpi=300)
        

	#------------------three plots side-by-side----------------

	histogram!(pl_b, L, normalize=:pdf, bins=100,
		color=color_arr[i],
		alpha=0.3,
		lw = 0,
		la=0.001,
		xlabel=L"L(t)",
		ylabel=ylabel,
	   	guidefontsize=guidefontsize,
		tickfontsize=tickfontsize,
		legendfontsize=legendfontsize,
		#label=L"\rho_{%$(i)}=%$(ρ)",
		label=L"b=%$(b)",
		legend=:topleft,
		framestyle=:box,
		ylims=[0,1.8e-2],
		yticks=yticks,
		xticks = xticks_forall[i],
		dpi=300,
		margin=5Plots.mm)
    

    plot!(pl_b, L_range, gauss.(L_range,μ,σ),
        color=color_arr[i],
        linewidth=2,
       	guidefontsize=guidefontsize,
		tickfontsize=tickfontsize,
		legendfontsize=legendfontsize,
        label=L"\mathcal{N}(S_{SK} \cdot t,\sqrt{\Lambda \cdot t})",
        #label=nothing,
        legend=:best,
        framestyle=:box,
        dpi=300)

	entropy_arrow_label = i == 1 ? L"S_{SK} \cdot t" : ""
	lyapunov_arrow_label = i == 1 ? L"\sqrt{\Lambda \cdot t}" : ""

    plot!(pl_b,[μ,μ],[0,gauss(μ,μ,σ)],lw=3,lc=logocolors.blue,arrow=:arrow,label=entropy_arrow_label)

	plot!(pl_b,[μ,μ+σ],[gauss(μ+σ,μ,σ),gauss(μ+σ,μ,σ)],lw=3,lc=logocolors.red,arrow=:arrow,label=lyapunov_arrow_label) 
     
	
    push!(individual_plots,pl_b)
    #savefig(pl_b, "walk_length_distribution_b=$(b)_Nsteps=$(N_steps)_ensemble=$(ensemble).svg")

end

annotate!(pl,8200, 0.009, text(L"S_{SK} = "*"$(round(entropies[1];digits=2))" * "\n" * L"\Lambda = "*"$(round(lyaps[1];digits=2))", :center, 10))
annotate!(pl,2600, 0.005, text(L"S_{SK} = "*"$(round(entropies[2];digits=2))" * "\n" * L"\Lambda = "*"$(round(lyaps[2];digits=2))", :center, 10))
annotate!(pl,6300, 0.011, text(L"S_{SK} = "*"$(round(entropies[3];digits=2))" * "\n" * L"\Lambda = "*"$(round(lyaps[3];digits=2))", :center, 10))

plot!(pl, xlim=(1900,9000))
savefig(pl, "walk_length_distribution_Nsteps=$(N_steps)_ensemble=$(ensemble).pdf")


l = @layout [a{0.33w} b{0.33w} c{0.33w}]

plot_all = plot(individual_plots...;layout=l,size=(1200,500),margin=7Plots.mm)

savefig(plot_all,"walklength_distributions_bs.pdf")






