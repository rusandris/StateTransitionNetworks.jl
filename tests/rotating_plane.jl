using StateTransitionNetworks
using DynamicalSystems
using ChaosTools
using Plots
using LinearAlgebra
using StatsBase
using LaTeXStrings
using DelimitedFiles
using Graphs
include("adding_stn_functions.jl")


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

function translate_to_origin!(traj::Matrix)
	for i in size(traj)[2]
		traj[:,i] = traj[:,i] .- mean(traj[:,i])
	end
end


#--------------------------------CHAOSTOOLS-----------------------------

_initialize_output(::S, ::Int) where {S} = eltype(S)[]
_initialize_output(u::S, i::SVector{N, Int}) where {N, S} = typeof(u[i])[]
function _initialize_output(u, i)
    error("The variable index when producing the PSOS must be Int or SVector{Int}")
end

function my_poincaresos(A::Dataset, plane; direction = -1, warning = true, idxs = 1:size(A, 2))
    i = typeof(idxs) <: Int ? idxs : SVector{length(idxs), Int}(idxs...)
    planecrossing = PlaneCrossing(plane, direction > 0)
    data = my_poincaresos(A, planecrossing, i)
    warning && length(data) == 0 && @warn PSOS_ERROR
    return Dataset(data)
end
function my_poincaresos(A::Dataset, planecrossing::PlaneCrossing, j)
    i, L = 1, length(A)
    data = _initialize_output(A[1], j)
    # Check if initial condition is already on the plane
    planecrossing(A[i]) == 0 && push!(data, A[i][j])
    i += 1
    side = planecrossing(A[i])

    while i ≤ L # We always check point i vs point i-1
        while side < 0 # bring trajectory infront of hyperplane
            i == L && break
            i += 1
            side = planecrossing(A[i])
        end
        
        # It is now guaranteed that A crosses hyperplane between i-1 and i
        ucross = interpolate_crossing(A[i-1], A[i], planecrossing)
        push!(data, ucross[j])
        
        while side ≥ 0 # iterate until behind the hyperplane
            i == L && break
            i += 1
            side = planecrossing(A[i])
        end
        i == L && break
        # It is now guaranteed that A crosses hyperplane between i-1 and i
        ucross = interpolate_crossing(A[i-1], A[i], planecrossing)
        push!(data, ucross[j])
    end
    return data
end

function interpolate_crossing(A, B, pc::PlaneCrossing{<:AbstractVector})
    # https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    t = LinearAlgebra.dot(pc.n, (pc.p₀ .- A))/LinearAlgebra.dot((B .- A), pc.n)
    return A .+ (B .- A) .* t
end

function interpolate_crossing(A, B, pc::PlaneCrossing{<:Tuple})
    # https://en.wikipedia.org/wiki/Linear_interpolation
    y₀ = A[pc.plane[1]]; y₁ = B[pc.plane[1]]; y = pc.plane[2]
    t = (y - y₀) / (y₁ - y₀) # linear interpolation with t₀ = 0, t₁ = 1
    return A .+ (B .- A) .* t
end

function rotplane_measures(traj;grid_size,plane0,θ_min,θ_max,Δθ,direction=+1)
	

	translate_to_origin!(traj)
	plane = deepcopy(plane0)
	rotate_plane!(plane,θ_min)
	θ = θ_min
	
	
	entropies = []
	lyapunovs = []
	PSOS_points_numbers = []
	average_degrees = []
	angles = []
	qweights = []
	
	while θ < θ_max

		if (direction ==:both) 
			psection =	my_poincaresos(traj, plane; save_idxs=[1,2,3],warning=true,direction=+1);
	
		else
			psection = DynamicalSystemsBase.poincaresos(traj, plane; save_idxs=[1,2,3],warning=true,direction=direction);
		end
		
		nrPSOS_points = length(psection)

		if nrPSOS_points <= 1
			@warn "Change the plane parameters!"
			
		end

		transform_to_plane!(psection,θ)
		traj_grid, vertex_names = timeseries_to_grid(psection[:,1:2],grid_size) 
		
		stn,retcode = create_stn(traj_grid,vertex_names,make_ergodic=true,verbose=true)
		if retcode ==:Success
			P = prob_matrix(stn)
			entropy, lyap = network_measures(P)
			
			push!(entropies,entropy)
			push!(lyapunovs,lyap)
			push!(PSOS_points_numbers,nrPSOS_points)
			push!(average_degrees,mean(outdegree(stn)))
			push!(qweights,weight_matrix(stn).nzval...)
		else
			push!(entropies,NaN)
			push!(lyapunovs,NaN)
			push!(PSOS_points_numbers,NaN)
			push!(average_degrees,NaN)
		end
		
		rotate_plane!(plane,Δθ)
		
		
		push!(angles,θ)
		
		
		θ += Δθ
		
	end
	return angles,entropies,lyapunovs,PSOS_points_numbers,average_degrees,qweights
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

