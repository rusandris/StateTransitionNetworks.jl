using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs


guidefontsize=20
tickfontsize=15
legendfontsize=15


include("rotating_plane.jl")

T = 5000
Ttr = 500
b_vals = [0.42,0.372,0.28,0.34]
ds = Systems.roessler();
grid_size = 20
direction = +1

plane0 = [0.0,0.0,1.0,0.0]
plots = []

for (i,b) in enumerate(b_vals)

	set_parameter!(ds,2,b)
	timeseries,  = trajectory(ds, T; Δt=0.01, Ttr=Ttr);
	
	@show b
	
	angles,entropies,lyaps,PSOS_points_numbers,average_degrees = rotplane_measures_v3(timeseries;
		grid_size=grid_size,
		plane0=plane0,
		θ_min=0,
		θ_max=π,
		Δθ=0.001,
		direction=direction,return_angles =true)
		
	#@show length(angles)
	#@show length(lyaps)
	#@show length(entropies)
	
	xlabel = i == 3 ? L"\theta" : ""
	
	pl_b = plot(angles,lyaps,
	dpi=300,
	c=:red,
	ylims = (0,1.3),
	yticks = [0.0,0.5,1.0],
	guidefontsize=guidefontsize,
	legendfontsize=legendfontsize,
	tickfontsize=tickfontsize,
	label=L"\Lambda",
	xlabel=xlabel,
	ylabel=L"~S_{KS},~\Lambda",
	framestyle=:box)


	plot!(pl_b,angles,entropies,
		tickfontsize=12,
		dpi=300,
		c=:blue,
		label=L"S_{KS}",
		framestyle=:box)
		
	i != 4 && push!(plots,pl_b)
		
	#=
	plot!(angles,PSOS_points_numbers/T,
		tickfontsize=12,
		dpi=300,
		c=:orange,
		guidefontsize=15,
		legendfontsize=12,
		label=L"n_{PSOS}/T",
		xlabel=L"grid ~ size",
		ylabel=L"n_{PSOS}/T,~S_{KS},~\Lambda",
		framestyle=:box)


	plot!(angles,average_degrees,
		tickfontsize=12,
		dpi=300,
		c=:green,
		guidefontsize=15,
		legendfontsize=12,
		label=L"<k>",
		xlabel=L"\theta",
		ylabel=L"n_{PSOS}/T,~S_{KS},~\Lambda,~<k>",
		framestyle=:box)
	=#
	    savefig(pl_b, "rotating_plane_roessler_b_$(b)"*"_dir$direction"*"T_$T"*".svg")
	    
	
end



l = @layout [a{0.33h}; b{0.33h}; c{0.33h}]

plot_all = plot(plots...;layout=l,size=(1000,1000),margin=5Plots.mm)

savefig(plot_all,"rotplane_measures_roessler.pdf")



