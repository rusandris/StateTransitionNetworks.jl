using DelimitedFiles
using Plots,LaTeXStrings

guidefontsize=20
tickfontsize=15
legendfontsize=15



#read data

data = readdlm("data_reproduced_netmeasures_roessler_b_T_5000_Ttr_500_traj_ens_1_rwens_100_nsteps_100_grid_20standardpsos.txt")

lyap_expdata = readdlm("lyapunov_exponent_roessler_b_T_5000_Ttr_500.txt")

od_data = readdlm("orbit_diagram_roessler_b_saved_z_T_5000_Ttr_500.txt")

parameter_values = data[:,1]
analytic_entropies = data[:,2]
analytic_lyapunovs = data[:,3]
num_entropies = data[:,4]
num_lyapunovs = data[:,5]


#---------------------orbit diagrams-------------------------------


special_bs = [0.42,0.369,0.34,0.28]
colors = [colorant"orange",colorant"red",colorant"blue",colorant"green"]

od_plot = plot(guidefontsize=guidefontsize,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	framestyle=:box,
	ylims=(-10,-2),
	ylabel=L"x",
	dpi=300,
	size=(800,600),
	reuse=false,
	legend=false)
	
for (i, p) in enumerate(parameter_values)
    od_line = od_data[i,od_data[i,:] .!= ""]
    
    plot!(od_plot,fill(p,length(od_line)),od_line,
    st=:scatter,
    ms=0.5,
    markerstrokewidth=0.00001,
    mc =:gray10)	
end
plot!(special_bs,fill(-2.5,4),st=:scatter,mc=colors,ms=8,markerstrokewidth=0.001)
#plot!([special_bs[2],special_bs[2]],[-2,-10.0],lc=:red,lw=0.3)
#plot!([special_bs[4],special_bs[4]],[-2,-10.0],lc=:green,lw=0.3)


plot!(title="Rössler system",titlefontsize=20,reuse=false)

#savefig(od_plot,"orbit_diagram_rössler.png")
    
#---------------------lyapunov exponent-------------------------------

lyap_exponents = lyap_expdata[:,2]

pl_lyap_exp = plot(parameter_values,lyap_exponents,
	lc=:gray10,
	ylabel=L"\lambda",
	ylims = (0,0.1),
	guidefontsize=guidefontsize,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	framestyle=:box,
	legend=false,
	dpi=300,
	lw = 2,
	yticks=[0.0,0.05,0.1],
	size=(800,300))

#savefig(pl_lyap_exp,"lyap_exponents_rössler.png")



#---------------------network measures-------------------------------



pl_lyap = plot(parameter_values,analytic_lyapunovs,
	lc=:gray10,
	lw=2,
	ylabel=L"\Lambda",
	guidefontsize=guidefontsize,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	framestyle=:box,
	label= "analytic",
	dpi=300,
	yticks=[0.0,0.8,1.6],
	size=(800,300))
	
plot!(pl_lyap,parameter_values,num_lyapunovs,
	mc=:gray50,
	ms=4,
	ma=0.4,
	markerstrokewidth = 0.001,
	st=:scatter,
	label= "random walk")	
	

#savefig(pl_lyap,"lyapunov_network_measure_rössler.svg")

pl_entr = plot(parameter_values,analytic_entropies,
	lc=:gray10,
	lw=2,
	xlabel=L"b",
	ylabel=L"S_{SK}",
	ylims=(0.0,1.2),
	guidefontsize=guidefontsize,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	framestyle=:box,
	label= "analytic",
	yticks=[0.0,0.6,1.2],
	dpi=300,
	size=(800,300))
	
plot!(pl_entr,parameter_values,num_entropies,
	mc=:gray50,
	ms=4,
	ma=0.4,
	markerstrokewidth = 0.001,
	st=:scatter,
	label= "random_walk",reuse=false)

#savefig(pl_entr,"entropies_network_measure_rössler.svg")

l = @layout [a{0.4h}; b{0.2h}; c{0.2h}; b{0.2h}]


plot_all = plot(od_plot,pl_lyap_exp,pl_lyap,pl_entr,size=(1000,1000),layout=l)

savefig(plot_all,"roesler_results.png")

