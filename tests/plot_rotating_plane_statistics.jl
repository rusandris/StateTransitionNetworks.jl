using Plots

PSOS_statistics = readdlm("network_measures_nrpsos_20_T_1000_Ttr_500.txt")
lyap_statistics = readdlm("network_measures_lyapunovs_20_T_1000_Ttr_500.txt")
entr_statistics = readdlm("network_measures_entropies_20_T_1000_Ttr_500.txt")
rho_values = PSOS_statistics[:,1]
angles = 0:0.001:π/2
T = 1000
		
for (i,rho) in enumerate(rho_values)
	#=
	hvert = histogram(PSOS_statistics[i,2:end] ./ T,
		bins=30,
		normalize=:probability,
		tickfontsize=6,
		title = "ρ = $rho",
		dpi=300,
		c=:purple,
		guidefontsize=15,
		legend=false,
		xlabel=L"n_{PSOS}/T",
		ylabel=L"P(n_{PSOS}/T)")
		
	push!(hists,hvert)
	=#
	angle_plot = plot(angles,PSOS_statistics[i,2:end] ./ T,
		tickfontsize=12,
		title = "ρ = $rho",
		dpi=300,
		c=:green,
		guidefontsize=15,
		legendfontsize=15,
		label=L"n_{PSOS}/T",
		xlabel=L"\theta",
		ylabel=L"n_{PSOS}/T,~S_{KS},~\Lambda",
		framestyle=:box)
		
	plot!(angles,lyap_statistics[i,2:end],
		tickfontsize=12,
		title = "ρ = $rho",
		dpi=300,
		c=:red,
		lw=2,
		guidefontsize=15,
		label=L"\Lambda",
		xlabel=L"\theta",
		ylabel=L"n_{PSOS}/T,~S_{KS},~\Lambda",
		framestyle=:box)
		
				
	plot!(angles,entr_statistics[i,2:end],
		tickfontsize=12,
		title = "ρ = $rho",
		dpi=300,
		c=:blue,
		lw=2,
		guidefontsize=15,
		label=L"S_{KS}",
		xlabel=L"\theta",
		ylabel=L"n_{PSOS}/T,~S_{KS},~\Lambda",
		framestyle=:box)
		
		savefig("rotating_plane_statistics_rho$rho"*"0_pi2"*".pdf")
end




