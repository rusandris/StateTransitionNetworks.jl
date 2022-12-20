using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using GraphPlot
using Random
using Statistics

using LinearAlgebra

# Lyapunov measure for a given ρ value
Δt = 0.001;
plane = (1,15.0);
grid = 20;
ρ=180.867;
u0 = rand(Float64,3).*50 .-25;
T = 3000;
system = Systems.lorenz(u0; ρ=ρ);
timeseries = trajectory(system, T; Δt=Δt, Ttr=500);

psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
d_traj, v_names = timeseries_to_grid(psection, grid);
stn, ret_code = create_stn(d_traj, v_names);
ret_code

P = prob_matrix(stn);
Q = weight_matrix(stn);

λ, v = eigen(Matrix(P));
λ

λ2 = (1 .-λ)

λ2 = λ2.^-1

#This needs better solution
λ2[abs.(λ2).>1e10] .= 0.0 + 0.0im;

Ω = Diagonal(λ2)


L = Matrix(-log.(P));
replace!(L, Inf=>0.0);
L = P.*L;

covariance = (v[:,end]./sum(v[:,end]))'*L*v*Ω*inv(v)*L*ones(length(λ))


avg = sinai_kolmogorov_entropy(Q,P)

sq_avg = sum(Q[Q .!=0] .* (log.(P[P .!=0])).^2)

sq_avg - avg^2 + 2*covariance

# Random walk Lyapunov
ensemble = 100;
N_steps = 10000;
l = 0;
for i in 1:20
   @show l
   S, L = network_measures(stn, ensemble, N_steps)
   l += L
end
l /= 20

# Analytic Lyapunov measure over an interval ======================================================================================

rho = 180:0.01:182;
T = 3000;
ensemble = 100;
data = [];
data2 = [];
data3 = [];
for ρ in rho
   @show ρ
   while true # V matrix should not be singular
      while true # STN must be 'healthy'
         u0 = rand(Float64,3).*50 .-25;
         system = Systems.lorenz(u0; ρ=ρ);
         timeseries = trajectory(system, T; Δt=Δt, Ttr=500);
         psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
         d_traj, v_names = timeseries_to_grid(psection, grid);
         stn, ret_code = create_stn(d_traj, v_names);
         if ret_code == :Success
            break
         end
      end

      P = prob_matrix(stn);
      Q = weight_matrix(stn);

      λ, v = eigen(Matrix(P));
      if det(v) != 0.0+0.0im
         break
      end
   end

   λ2 = (1 .-λ)

   λ2 = λ2.^-1

   # This needs better solution
   λ2[abs.(λ2).>1e10] .= 0.0 + 0.0im;

   Ω = Diagonal(λ2)

   L = Matrix(-log.(P));
   replace!(L, Inf=>0.0);
   L = P.*L;

   covariance = (v[:,end]./sum(v[:,end]))'*L*v*Ω*inv(v)*L*ones(length(λ))
   avg = sinai_kolmogorov_entropy(Q,P)

   sq_avg = sum(Q[Q .!=0] .* (log.(P[P .!=0])).^2)

   lyapunov = real(sq_avg - avg^2 + 2*covariance)
   push!(data, lyapunov)

   S, L = network_measures(stn, ensemble, N_steps)
   push!(data2, L)
   push!(data3, S)
   
end

plot(rho, data, label = "Analytic", xlabel = "ρ", ylabel = "Λ")

plot!(rho, data2, label = "Random walk")
plot!(legend = :bottomleft)
