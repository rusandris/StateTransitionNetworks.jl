using DynamicalSystems
using StateTransitionNetworks
using Plots
using LaTeXStrings


guidefontsize=20
tickfontsize=15
legendfontsize=15

ds = Systems.roessler()
special_bs = [0.42,0.368,0.34,0.28]
ensembles = 10 .^[2,3,4]
ensemble = ensembles[3]
N_max = 10000
plane = (2,0.0)
T = 5000
Ttr = 500
steps = 1:N_max
colors = [:orange,:red,:blue,:green]
grid_size = 20

#-----------------------for 3 ensemble sizes------------------------

alphas = [0.2,0.6,1]
pl_lyap = plot()
pl_entr = plot()
b = special_bs[1]

set_parameter!(ds,2,b)
for (i,ens) in enumerate(ensembles)
	

	@show ens
	traj = trajectory(ds, T; Δt=0.01, Ttr=Ttr);
	stn,retcode = create_stn(traj,grid_size,plane,[1,3];direction=1)


	if retcode ==:Success
		entr,lyaps = measure_convergence(stn,ens,N_max)
		P = prob_matrix(stn)
		analytic_entr,analytic_lyap = network_measures(P)
		
			
		plot!(pl_lyap,lyaps,
		lw=2,
		ylims = (0.2,0.5),
		yticks = [0.2,0.35,0.5],
		la = alphas[i],
		lc =:gray10,
		xlabel=L"t",
		ylabel=L"\Lambda",
		guidefontsize=guidefontsize,
		tickfontsize=tickfontsize,
		legendfontsize=legendfontsize,
		label=L"N = "*"$ens",
		framestyle=:box,
		legend=false,
		xticks = [0,Int(N_max/2),N_max],
		title = L"b = " *"$b",
		dpi=300)
		
		plot!(pl_lyap,[steps[1],steps[end]],[analytic_lyap,analytic_lyap],lw = 1, ls=:dash,lc=:gray10,label="")
		
		
		
					
		plot!(pl_entr,entr,
		lw=2,
		ylims = (0.8,0.85),
		yticks = [0.8,0.825,0.85],
		la = alphas[i],
		lc = :gray10,
		xlabel=L"t",
		ylabel=L"S_{SK}",
		guidefontsize=guidefontsize,
		tickfontsize=tickfontsize,
		legendfontsize=legendfontsize,
		label=L"N = "*"$ens",
		legend=:false,
		title = L"b = " *"$b",
		framestyle=:box,
		xticks = [0,Int(N_max/2),N_max],
		dpi=300)
		
		plot!(pl_entr,[steps[1],steps[end]],[analytic_entr,analytic_entr],lw = 1, ls=:dash,lc=:gray10,label="")
		
	end

end

empty_plot = plot((1:3)',legend=:outertopleft,legendfontsize=legendfontsize,color = [:gray10,:gray10,:gray10],la=alphas',labels=[L"N=10^2" L"N=10^3" L"N=10^4"],framestyle = :none)

pl_ensembles = plot(pl_entr,pl_lyap,empty_plot,layout=(1,3),size=(1000,300),margin=9Plots.mm)

savefig(pl_ensembles,"roessler_convergence_b_$b"*"_N$N_max"*"_T$T"*"_Ttr$Ttr.pdf")





#-------------------------for one ensemble, all params-------------------


pl_lyap = plot()
pl_entr = plot()

for (i,b) in enumerate(special_bs)
	set_parameter!(ds,2,b)
	@show b
	traj = trajectory(ds, T; Δt=0.01, Ttr=Ttr);
	stn,retcode = create_stn(traj,grid_size,plane,[1,3],direction=1)
	
	if retcode ==:Success
		entr,lyaps = measure_convergence(stn,ensemble,N_max)
		P = prob_matrix(stn)
		analytic_entr,analytic_lyap = network_measures(P)
		
		plot!(pl_lyap,lyaps,
		lw=1,
		lc = colors[i],
		xlabel=L"t",
		ylabel=L"\Lambda",
		guidefontsize=guidefontsize,
		tickfontsize=tickfontsize,
		legendfontsize=legendfontsize,
		title = L"N = 10^4",
		framestyle=:box,
		ylims=(0.0,1.8),
		yticks = [0.0,0.9,1.8],
		xticks = [0,Int(N_max/2),N_max],
		legend=false,
		dpi=300)
		
		plot!(pl_lyap,[steps[1],steps[end]],[analytic_lyap,analytic_lyap],lw = 1, ls=:dash,lc=colors[i],ma=0.4,label="")
		
		plot!(pl_entr,entr,
		lw=1,
		lc = colors[i],
		xlabel=L"t",
		ylabel=L"S_{SK}",
		guidefontsize=guidefontsize,
		tickfontsize=tickfontsize,
		legendfontsize=legendfontsize,
		label=L"b = "*"$b",
		title = L"N = 10^4",
		framestyle=:box,
		ylims=(0.0,1.0),
		yticks = [0.0,0.5,1.0],
		xticks = [0,Int(N_max/2),N_max],
		legend=false,
		dpi=300)
		
		plot!(pl_entr,[steps[1],steps[end]],[analytic_entr,analytic_entr],lw = 1, ls=:dash,lc=colors[i],ma=0.4,label="")
		
		
		
	end
end

empty_plot = plot((1:4)',legend=:outertopleft,legendfontsize=legendfontsize,color = [:orange :red :blue :green] ,labels=[L"b=0.42" L"b=0.368" L"b=0.34" L"b=0.28"],framestyle = :none)
pl_bs = plot(pl_entr,pl_lyap,empty_plot,layout=(1,3),size=(1000,300),margin=9Plots.mm)

savefig(pl_bs,"roessler_convergence_ens$ensemble"*"_N$N_max"*"_T$T"*"_Ttr$Ttr.pdf")
plot_all = plot(pl_ensembles,pl_bs,layout=(2,1),size=(1000,600))
savefig(plot_all,"roessler_convergence_T$T"*"_Ttr$Ttr.pdf")



