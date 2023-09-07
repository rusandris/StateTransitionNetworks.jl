using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs


guidefontsize=20
tickfontsize=15
legendfontsize=15


include("../rotating_plane.jl")

T = 5000
Ttr = 500
rho_vals = [180.1, 180.7,180.78]
ds = PredefinedDynamicalSystems.lorenz();
grid_size = 20
direction = +1

plane0 = [0.0,0.0,1.0,0.0]
plots = []

for (i,rho) in enumerate(rho_vals)

	set_parameter!(ds,2,rho)
	timeseries,  = trajectory(ds, T; Δt=0.01, Ttr=Ttr);
	
	@show rho
	
	angles,entropies,lyaps,PSOS_points_numbers,average_degrees = rotplane_measures(timeseries;
		grid_size=grid_size,
		plane0=plane0,
		θ_min=0,
		θ_max=π,
		Δθ=0.01,
		direction=direction)
		
	writedlm("data/data_rotating_plane_lorenz_rho_$(rho)"*"_dir$direction"*"T_$T"*"grid_$grid_size"*".txt",hcat(angles,entropies,lyaps,PSOS_points_numbers,average_degrees))
		
	#@show length(angles)
	#@show length(lyaps)
	#@show length(entropies)
	
	xlabel = i == 3 ? L"\theta" : ""
	
	pl_b = plot(angles,lyaps,
	dpi=300,
	c=:red,
	guidefontsize=guidefontsize,
	legendfontsize=legendfontsize,
	tickfontsize=tickfontsize,
	label=L"\Lambda",
	xlabel=xlabel,
	ylabel=L"~S_{KS},~\Lambda",
	framestyle=:box)


	plot!(pl_b,angles,entropies,
		dpi=300,
		c=:blue,
		label=L"S_{KS}",
		framestyle=:box)
		
	i != 4 && push!(plots,pl_b)
		

	
	pl_other = plot(angles,PSOS_points_numbers/T,
		dpi=300,
		c=:orange,
		guidefontsize=guidefontsize,
		legendfontsize=legendfontsize,
		tickfontsize=tickfontsize,
		label=L"n_{PSOS}/T",
		ylabel=L"n_{PSOS}/T,~S_{KS},~\Lambda",
		framestyle=:box)


	plot!(pl_other,angles,average_degrees,
		dpi=300,
		c=:green,
		guidefontsize=guidefontsize,
		legendfontsize=legendfontsize,
		tickfontsize=tickfontsize,
		label=L"<k>",
		xlabel=L"\theta",
		ylabel=L"n_{PSOS}/T,~S_{KS},~\Lambda,~<k>",
		framestyle=:box)
		
		l = @layout [a{0.5h}; b{0.5h}]
		plot_all_measures = plot(pl_b,pl_other;layout=l,size=(1000,500),margin=5Plots.mm)
	
	    savefig(plot_all_measures, "figs/rotating_plane_lorenz_rho_$(rho)"*"_dir$direction"*"T_$T"*".pdf")
	    
	
end



l = @layout [a{0.33h}; b{0.33h}; c{0.33h}]

plot_all = plot(plots...;layout=l,size=(1000,1000),margin=5Plots.mm)

savefig(plot_all,"figs/rotplane_measures_lorenz.pdf")



