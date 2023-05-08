using StateTransitionNetworks
using DynamicalSystems
using ChaosTools
using Plots
using LinearAlgebra
using StatsBase
using LaTeXStrings
using DelimitedFiles

#rotation matrix - counterclockwise around y
#theta measured between normal and z axis
Ry(θ) = [cos(θ) 0 sin(θ);0 1 0; -sin(θ) 0 cos(θ)]

#poincare plane given as plane = [a1,a2,a3,b]
#where a1*x1+a2*x2+a3*x3 = b
function rotate_plane!(plane,θ)
	a = plane[1:end-1]
	b = plane[end]
	
	a = Ry(θ)*a
	
	plane[:] = [a...,b]
end

#transform to in-plane coordinates
function transform_to_plane!(psection,θ)
	
	for i in 1:length(psection)
		psection[i] = Ry(-θ)*psection[i]
	end

end

function translate_to_origin!(traj)
	dim = length(traj[1])
	means = [mean(traj[:,i]) for i in 1:dim]
	
	for s in 1:length(traj.data)
		traj.data[s] = traj.data[s] - means
	end
end

function parameter_sweep(rho_values;grid_size,T,Ttr,Δt,θ_min,θ_max,Δθ,save = false)
	ds = Systems.lorenz()


	plane0 = [0.0,0.0,1.0,0.0]
	
	nrPSOS_points_plot = plot(ms=10,
	xlabel = L"$\rho$",
	ylabel = "number of PSOS points",
	xguidefontsize=20,
	yguidefontsize=20,
	tickfontsize=10,
	alpha = 0.5,
	legend=false,
	dpi = 300)

	lyap_plot = plot(ms=10,
	xlabel = L"$\rho$",
	ylabel = L"$\Lambda$",
	xguidefontsize=20,
	yguidefontsize=20,
	tickfontsize=10,
	alpha = 0.5,
	legend=false,
	dpi = 300)
	
	
	entr_plot = plot(ms=10,
	xlabel = L"$\rho$",
	ylabel = L"$S_{SK}$",
	xguidefontsize=20,
	yguidefontsize=20,
	tickfontsize=10,
	markerstrokewidth=0.0001,
	alpha = 0.5,
	legend=false,
	dpi = 300)
	
	if save
		lyap_container = []
		entr_container = []
		nrPSOS_container = []
			
		entropy_averages = []
		lyapunov_averages = []
		nrPSOS_points_averages = []
	end

	for rho in rho_values
		@show rho
		set_parameter!(ds,2,rho)
		traj = trajectory(ds,T;Ttr = Ttr,Δt = 1e-3)
		translate_to_origin!(traj)
		
		entropies,lyapunovs,nrPSOS_points = rotplane_measures(traj;grid_size,plane0,θ_min,θ_max,Δθ)
		push!(entropy_averages,mean(entropies))
		push!(lyapunov_averages,mean(lyapunovs))
		push!(nrPSOS_points_averages,mean(nrPSOS_points))
		
		if save
			push!(lyap_container,lyapunovs)
			push!(entr_container,entropies)
			push!(nrPSOS_container,nrPSOS_points)
		end

		plot!(lyap_plot,fill(rho,length(lyapunovs)),lyapunovs,st=:scatter,mc=:red,alpha=0.3,markerstrokewidth=0.0001)
		plot!(entr_plot,fill(rho,length(entropies)),entropies,st=:scatter,mc=:blue,alpha=0.3,markerstrokewidth=0.0001)
		plot!(nrPSOS_points_plot,fill(rho,length(entropies)),nrPSOS_points,st=:scatter,mc=:blue,alpha=0.3,markerstrokewidth=0.0001)
	end
	plot!(lyap_plot,rho_values,lyapunov_averages,lw=1,lc=:gray10,legend=false)
	plot!(entr_plot,rho_values,entropy_averages,lw=1,lc=:gray10,legend=false)
	plot!(nrPSOS_points_plot,rho_values,nrPSOS_points_averages,lw=1,lc=:gray10,legend=false)
	
	if save 
		writedlm("network_measures_lyapunovs_$grid_size"*"_T_$T"*"_Ttr_$Ttr"*".txt",hcat(rho_values,hcat(lyap_container...)'))
		writedlm("network_measures_entropies_$grid_size"*"_T_$T"*"_Ttr_$Ttr"*".txt",hcat(rho_values,hcat(entr_container...)'))
		writedlm("network_measures_nrpsos_$grid_size"*"_T_$T"*"_Ttr_$Ttr"*".txt",hcat(rho_values,hcat(nrPSOS_container...)'))
		savefig(lyap_plot,"network_measures_lyapunovs_$grid_size"*"_T_$T"*"_Ttr_$Ttr"*".png")
		savefig(entr_plot,"network_measures_entropies_$grid_size"*"_T_$T"*"_Ttr_$Ttr"*".png")
	else
		return lyap_plot,entr_plot
	end
end


function rotplane_measures(traj;grid_size,plane0,θ_min,θ_max,Δθ,return_angles=false)
	
	entropies = []
	lyapunovs = []
	PSOS_points_numbers = []
	
	plane = deepcopy(plane0)
	rotate_plane!(plane,θ_min)
	θ = θ_min
	
	if return_angles
		angles = []
	end
	
	while θ < θ_max

		psection = ChaosTools.poincaresos(traj, plane; idxs=[1,2,3],warning=true,direction=+1);
		nrPSOS_points = length(psection)
		@show nrPSOS_points
		if nrPSOS_points <= 1
			@warn "Change the plane parameters!"
			continue
		end
		
		transform_to_plane!(psection,θ)
		traj_grid, vertex_names = timeseries_to_grid(psection[:,1:2],grid_size) 
		
		stn,retcode = create_stn(traj_grid,vertex_names)
		if retcode ==:Success
			P = prob_matrix(stn)
			entropy, lyap = network_measures(P)
			
			push!(entropies,entropy)
			push!(lyapunovs,lyap)
			push!(PSOS_points_numbers,nrPSOS_points)
		else
			push!(entropies,NaN)
			push!(lyapunovs,NaN)
		end
		
		rotate_plane!(plane,Δθ)
		
		if return_angles
			push!(angles,θ)
		end
		
		θ += Δθ
		
	end
	if return_angles
		return angles,entropies,lyapunovs
	else
		return entropies,lyapunovs,PSOS_points_numbers
	end
end
