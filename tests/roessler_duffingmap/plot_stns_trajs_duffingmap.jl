using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using LaTeXStrings
using StatsBase
using DelimitedFiles
using LinearAlgebra

T = 3000
Ttr = 1000
grid_size=20

guidefontsize=20
tickfontsize=15
legendfontsize=15


data = readdlm("data_netmeasures_duffingmap_b_T_30000_Ttr_1000_rwens_100_nsteps_100_grid_20.txt")

parameter_values = data[:,1]
analytic_lyapunovs = data[:,3]

od_data = readdlm("orbit_diagram_duffingmap_a_saved_x_T_30000_Ttr_1000.txt")

@show parameter_values[argmax(analytic_lyapunovs)]


special_as = [2.5054,2.50579,2.507,2.5084]
colors = [colorant"orange",colorant"red",colorant"blue",colorant"green"]


duffingmap(u,p,n) = SVector(u[2],-p[2]*u[1] +p[1]*u[2] - u[2]^3) 
p = [2.3,0.1]
u0 = [0.5,0.2]

ds = DiscreteDynamicalSystem(duffingmap,u0,p);


#=

od_plot = plot(guidefontsize=guidefontsize,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	framestyle=:box,
	ylabel=L"x",
	ylims=(0.8,1.0),
	dpi=300,
	size=(800,600),
	reuse=false,
	legend=false)
	
for (i, p) in enumerate(parameter_values)
    od_line = od_data[i,od_data[i,:] .!= ""]
    
    plot!(od_plot,fill(p,length(od_line)),od_line,
    st=:scatter,
    ms=0.5,
    markerstrokewidth=0.001,
    mc =:gray10)	
end


plot!(special_as,fill(0.98,4),st=:scatter,mc=colors,ms=8,markerstrokewidth=0.001)
plot!([special_as[2],special_as[2]],[0.8,1.0],lc=:red,lw=0.3)

=#


for a in special_as
	@show a 
	
	set_parameter!(ds,1,a)
	@show ds.p
	pred = predictability(ds;Ttr = 1000,λ_max = lyapunov(ds, 5000;Ttr=1000))
	@show pred
end


trajectories = []
stns = []

	
#----------------plotting an stn---------------------
for (i,a) in enumerate(special_as)
	println("Plotting stn for special a $a...")
	
	set_parameter!(ds,1,a)
	@show a	


	timeseries = trajectory(ds, T; Ttr=Ttr);

	ms = i == 3 ? 4 : 1.5
	yticks = i == 1 ? true : false
	ylabel = i == 1 ? L"y" : ""

	traj = plot(timeseries[:,1],timeseries[:,2],
	st=:scatter,
	mc=colors[i],
	ms = ms,
	markerstrokewidth=0.0001,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	guidefontsize=guidefontsize,
	framestyle=:box,
	yticks=yticks,
	xlabel=L"x",
	ylabel=ylabel,
	label=L"a = "*"$a",
	dpi=300,
	xlims=(0.6,1.5),
	ma=0.05,
	legend=:bottomleft,
	left_margin = 7Plots.mm,
	bottom_margin=7Plots.mm)
	
	#savefig(traj3d, "Tinkerbell_traj_b_$b.png")
	push!(trajectories,traj)
	
	
	d_traj, v_names = timeseries_to_grid(timeseries, grid_size);
	stn, retcode = create_stn(d_traj, v_names;make_ergodic=true,verbose=true);
	
	if retcode == :Success
		P = prob_matrix(stn)
		c = RGBA(colors[i])

		p_nodes = real(eigvecs(Matrix(P)')[:,end]) 
		p_nodes /= sum(p_nodes)


		plot_stn(stn;filename="stn_duffingmap_a_$a"*".pdf",
			nodesize=0.2,
			nodefillc=[RGBA(c.r,c.g,c.b,x) for x in p_nodes .^ (1/8)],
			linetype="curve",
			max_edgelinewidth=1.0,
			nodelabeldist=1.5, 
			nodelabelangleoffset=π/4,
			arrowlengthfrac=0.05)
			
		println("Done.")
		
	end

end

l = @layout [a{0.25w} b{0.25w} c{0.25w} d{0.25w}]

plot_all = plot(trajectories...,layout=l,size=(1200,400))

savefig(plot_all,"duffingmap_trajectories.png")
