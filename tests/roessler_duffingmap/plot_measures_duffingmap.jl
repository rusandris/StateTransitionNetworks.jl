using DelimitedFiles
using Plots,LaTeXStrings

guidefontsize=20
tickfontsize=15
legendfontsize=15


data = readdlm("data_netmeasures_duffingmap_b_T_30000_Ttr_1000_rwens_100_nsteps_10000_grid_20.txt")

parameter_values = data[:,1]
analytic_entropies = data[:,2]
analytic_lyapunovs = data[:,3]
num_entropies = data[:,4]
num_lyapunovs = data[:,5]


od_data = readdlm("orbit_diagram_duffingmap_a_saved_x_T_30000_Ttr_1000.txt")

lyap_expdata = readdlm("lyapunov_exponent_duffingmap_a_T_30000_Ttr_1000.txt")
lyap_exponents = lyap_expdata[:,2]


#---------------------orbit diagrams-------------------------------


special_as = [2.5054,2.50579,2.507,2.5084]
colors = [colorant"orange",colorant"red",colorant"blue",colorant"green"]

od_plot = plot(guidefontsize=guidefontsize,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	framestyle=:box,
	ylabel=L"x",
	ylims=(0.8,1.0),
	yticks=[0.8,0.9,1.0],
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
#plot!([special_as[2],special_as[2]],[0.8,1.0],lc=:red,lw=0.3)

plot!(title="Duffing map",titlefontsize=20,reuse=false)

#savefig(od_plot,"orbit_diagram_rössler.png")
    
#---------------------lyapunov exponent-------------------------------

lyap_exponents = lyap_expdata[:,2]

pl_lyap_exp = plot(parameter_values,lyap_exponents,
	lc=:gray10,
	ylabel=L"\lambda",
	guidefontsize=guidefontsize,
	tickfontsize=tickfontsize,
	legendfontsize=legendfontsize,
	framestyle=:box,
	legend=false,
	ylims=(-0.4,0.4),
	yticks=[-0.4,0.0,0.4],
	dpi=300,
	lw = 2,
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
	ylims=(0.0,1.0),
	yticks=[0.0,0.5,1.0],
	dpi=300,
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
	yticks=[0.0,0.5,1.0],
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

savefig(plot_all,"duffingmap_results.png")





