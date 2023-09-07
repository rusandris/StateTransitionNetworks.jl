using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using LaTeXStrings
using StatsBase
using DelimitedFiles
using LinearAlgebra

T = 5000
Ttr = 1000
plane = (2,0.0)
grid_size=20

guidefontsize=20
tickfontsize=15
legendfontsize=15


data = readdlm("data/data_netmeasures_roessler_b_T_5000_Ttr_500_traj_ens_rwens_100_nsteps_10000_grid_20standardpsos.txt")

parameter_values = data[:,1]
analytic_lyapunovs = data[:,3]

od_data = readdlm("data/orbit_diagram_roessler_b_saved_z_T_5000_Ttr_500.txt")


@show parameter_values[argmax(analytic_lyapunovs)]

special_bs = [0.42,0.368,0.34,0.28]
colors = [colorant"orange",colorant"red",colorant"blue",colorant"green"]

#=
od_plot = plot(guidefontsize=guidefontsize,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	framestyle=:box,
	ylabel=L"x",
	dpi=300,
	size=(800,600),
	reuse=false,
	legend=false)
	
for (i, p) in enumerate(parameter_values)
    od_line = od_data[i,od_data[i,:] .!= ""]
    
    plot!(od_plot,fill(p,length(od_line)),od_line,
    st=:scatter,
    ms=0.1,
    markerstokewitdh=0.001,
    mc =:gray50)	
end

plot!(special_bs,fill(-2.5,3),st=:scatter,mc=colors,ms=8,markerstrokewidth=0.001)
plot!([special_bs[2],special_bs[2]],[-9,-3])

=#

ds = PredefinedDynamicalSystems.roessler();


#predictability not implemented in v3

#=
for b in special_bs
	@show b 
	set_parameter!(ds,2,b)
	pred = predictability(ds;Ttr = 500,T_max = 1000,λ_max = lyapunov(ds, 5000;Ttr=1000))
	@show pred
end
=#

trajectories = []
psos_plots = []
	
#----------------plotting an stn---------------------
for (i,b) in enumerate(special_bs)
	println("Plotting stn for special b $b...")
	
	set_parameter!(ds,2,b)
	@show b	


	timeseries,  = trajectory(ds, T; Δt=0.01, Ttr=Ttr);
	
	yticks = i == 1 ? true : false
	ylabel = i == 1 ? L"y" : ""

	lw = i == 3 ? 2 : 0.1

	traj = plot(timeseries[:,1],timeseries[:,2],
	lc=colors[i],
	lw=lw,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	guidefontsize=guidefontsize,
	framestyle=:box,
	xlabel=L"x",
	ylabel=ylabel,
	yticks=yticks,
	label=L"b = "*"$b",
	aspect_ratio=1,
	legend=:topright,
	ylims=(-15,15),
	xlims=(-10,12),
	left_margin = 5Plots.mm,
	bottom_margin=7Plots.mm,
	dpi=300)
	
	plot!([-9,-1],[0.0,0.0],lw=3,ls=:dot,lc=:gray10,label="Poincaré plane")
	
	#savefig(traj3d, "Roessler_traj_b_$b.png")
	
	push!(trajectories,traj)
	
	psection = DynamicalSystemsBase.poincaresos(timeseries, plane; direction=+1, save_idxs=[1,3]);
	
	
	ms = i == 3 ? 4 : 1.5
	ylabel = i == 1 ? L"z" : ""
	
	psos_plot = plot(psection[:,1],psection[:,2],
	st=:scatter,
	mc=colors[i],
	ms=ms,
	ma=0.1,
	xlims = (-10,-2),
	ylims = (0.018,0.05),
	markerstrokewidth=0.000,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	guidefontsize=guidefontsize,
	legend=:topleft,
	framestyle=:box,
	xlabel=L"x",
	ylabel=ylabel,
	yticks=yticks,
	xticks=[-9,-6.5,-4],
	label=L"b = "*"$b",
	left_margin = 5Plots.mm,
	bottom_margin=7Plots.mm,
	dpi=300)
	
	push!(psos_plots,psos_plot)
	
	
	stn, retcode = create_stn(psection, grid_size;make_ergodic=true,verbose=true);
	
	if retcode == :Success
		P = prob_matrix(stn)
		c = RGBA(colors[i])

		p_nodes = real(eigvecs(Matrix(P)')[:,end]) 
		p_nodes /= sum(p_nodes)


		plot_stn(stn;filename="figs/stn_roessler_b_$b"*".pdf",
			nodesize=0.2,
			nodefillc=[RGBA(c.r,c.g,c.b,x) for x in p_nodes .^ (1/8)],
			linetype="curve",
			max_edgelinewidth=1.0,
			nodelabeldist=1.5,
			weight_exponent=1,
			nodelabelangleoffset=π/4,
			arrowlengthfrac=0.05)
			
		println("Done.")
		
	end

end


l = @layout [a{0.25w} b{0.25w} c{0.25w} d{0.25w};a{0.25w} b{0.25w} c{0.25w} d{0.25w}]

plot_all = plot(trajectories...,psos_plots...,layout=l,size=(1200,800))

savefig(plot_all,"figs/roessler_trajectories_psos.pdf")
