using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs

include("rotating_plane.jl")

T = 5000
rho_vals = [180.1, 180.7, 180.78]
ds = Systems.lorenz();
grid_size = 20

plane0 = [0.0,0.0,1.0,0.0]


for rho in rho_vals

	set_parameter!(ds,2,rho)
	timeseries = trajectory(ds, T; Δt=0.01, Ttr=500);
	
	@show rho
	
	angles,entropies,lyaps,PSOS_points_numbers,average_degrees = rotplane_measures(timeseries;
		grid_size=grid_size,
		plane0=plane0,
		θ_min=0,
		θ_max=π,
		Δθ=0.001,
		return_angles=true)
		
		
	pl_rho = plot(angles,lyaps,
	tickfontsize=12,
	dpi=300,
	c=:red,
	guidefontsize=15,
	legendfontsize=12,
	label=L"\Lambda",
	xlabel=L"grid ~ size",
	ylabel=L"n_{PSOS}/T,~S_{KS},~\Lambda",
	framestyle=:box)


	plot!(angles,entropies,
		tickfontsize=12,
		dpi=300,
		c=:blue,
		guidefontsize=15,
		legendfontsize=12,
		label=L"S_{KS}",
		xlabel=L"grid ~ size",
		ylabel=L"n_{PSOS}/T,~S_{KS},~\Lambda",
		framestyle=:box)
		
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
		
	    savefig(pl_rho, "rotating_plane_rho_$(rho).svg")
	
end
