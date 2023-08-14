using Plots,LaTeXStrings
using DelimitedFiles

data_filename = "data_rotating_plane_lorenz_rho_180.1_dir1T_5000grid_20"
rotplane_data = readdlm("data/"*data_filename*".txt")

rho = 180.1


guidefontsize=20
tickfontsize=15
legendfontsize=15


angles = rotplane_data[:,1]
entropies = rotplane_data[:,2]
lyaps = rotplane_data[:,3]
PSOS_points_numbers = rotplane_data[:,4]
average_degrees = rotplane_data[:,5]

#@show length(angles)
#@show length(lyaps)
#@show length(entropies)

pl_b = plot(angles,lyaps,
dpi=300,
c=:red,
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
label=L"\Lambda",
xlabel="",
ylabel=L"~S_{KS},~\Lambda",
framestyle=:box)


plot!(pl_b,angles,entropies,
	dpi=300,
	c=:blue,
	label=L"S_{KS}",
	ylims = (0.0,2.0),
	yticks = [0.0,1.0,2.0],
	framestyle=:box)
	

pl_other = plot(angles,PSOS_points_numbers/T,
	dpi=300,
	c=:orange,
	guidefontsize=guidefontsize,
	legendfontsize=legendfontsize,
	tickfontsize=tickfontsize,
	ylims = (1,7),
	yticks = [1,4,7],
	label=L"n_{PSOS}/T",
	ylabel=L"n_{PSOS}/T, \langle k \rangle",
	framestyle=:box)


plot!(pl_other,angles,average_degrees,
	dpi=300,
	c=:green,
	guidefontsize=guidefontsize,
	legendfontsize=legendfontsize,
	tickfontsize=tickfontsize,
	label=L"<k>",
	xlabel=L"\theta",
	framestyle=:box)
	
l = @layout [a{0.5h}; b{0.5h}]
plot_all_measures = plot(pl_b,pl_other;layout=l,size=(1000,500),margin=5Plots.mm)

savefig(plot_all_measures, "figs/"*data_filename*".pdf")
