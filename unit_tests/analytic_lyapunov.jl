using StateTransitionNetworks
using DynamicalSystems
using Test

ensemble = 10^3
T = 10
Δt = 0.001;
plane = (1,15.0);
grid_size = 20
L = []
rho = [180.1, 180.7, 180.78];

for i in eachindex(rho)
	@show ρ=rho[i]
	l = []

	while length(l) < ensemble
		u0 = rand(Float64,3).*50 .-25;
		system = PredefinedDynamicalSystems.lorenz(u0; ρ=ρ);
		timeseries,  = trajectory(system, T; Δt=Δt, Ttr=1000);
		psection = DynamicalSystemsBase.poincaresos(timeseries, plane; direction=+1, save_idxs=[2,3]);
	    d_traj, v_names = timeseries_to_grid(psection, grid_size);
	    stn, retcode = create_stn(d_traj, v_names; make_ergodic=true,verbose=false);
	    
	    
	    if retcode == :Success
		    P = prob_matrix(stn);
		    Q = weight_matrix(stn);
		    @test isnormalized(P)
		    if !any(isnan, P)
		        lyapunov = lyapunov_measure(P)
		        @test lyapunov[4] == :Success
		        #@test lyapunov >= 0
		        #@test (lyapunov > 0 || lyapunov ≈ 0)
		        push!(l, lyapunov[1])
		    end
		end
		
	end
	push!(L,l);
	
end
