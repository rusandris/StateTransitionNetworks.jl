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
ρ= 181.02; #180.9;
u0 = rand(Float64,3).*50 .-25;
# ρ = 180.91
# u0 = [19.208757974895896, 8.16404896961214, 22.003233673468387]
# # u0 = [-6.415632455803806, 22.288031139971572, -0.19849359973168035]
# # u0 = [-3.7172591949410076, 19.267269973658813, -10.673236379143875]
T = 3000;
system = Systems.lorenz(u0; ρ=ρ);
timeseries = trajectory(system, T; Δt=Δt, Ttr=500);

psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
d_traj, v_names = timeseries_to_grid(psection, grid);
stn, ret_code = create_stn(d_traj, v_names);
ret_code

P = prob_matrix(stn);
Q = weight_matrix(stn);

function analytic_lyapunov(P,Q)
   λ, v = eigen(Matrix(P))
   @show det(v)
   if isapprox(det(v),0) 
      return -1, -1, -1, :SingularityError
   end

   λt, vt = eigen(transpose(Matrix(P)))

   if real(λ[end]) ≈ 1
      λ[end] = 1.0 + 0.0im
   end

   if real(λt[end]) ≈ 1
      λt[end] = 1.0 + 0.0im
      x = (vt[:,end]./sum(vt[:,end]))'
   end

   Ω = Diagonal(1 ./ (1 .- λ))
   replace!(Ω, NaN+NaN*im=>0.0)

   L = Matrix(-log.(P))
   replace!(L, Inf=>0.0)
   L = P.*L

   v*inv(v)
   det(v)
   det(v)≈ 0.0+0.0im
   covariance = x*L*v*Ω*inv(v)*L*ones(length(λ))

   avg = sinai_kolmogorov_entropy(Q,P)
   sq_avg = sum(Q[Q .!=0] .* (log.(P[P .!=0])).^2)

   variance = sq_avg - avg^2
   lyapunov =  variance + 2*real(covariance)
   if imag(covariance) < 1.0e-3
      return lyapunov, variance, real(covariance), :Success
   else
      @show covariance
      return lyapunov, variance, real(covariance), :ComplexCovariancveWarning
   end
end

analytic_lyapunov(P,Q)
network_measures(stn, 1000, 10^4)[2]


# Analytic Lyapunov measure over an interval ======================================================================================

rho = 180:0.005:182;
T = 5000;
ensemble = 1000;
N_steps = 10000;
data0 = [];
data1 = [];
data2 = [];
data3 = [];
data4 = [];
for ρ in rho
      @show ρ
      while true # STN must be 'healthy'
         u0 = rand(Float64,3).*50 .-25;
         system = Systems.lorenz(u0; ρ=ρ);
         timeseries = trajectory(system, T; Δt=Δt, Ttr=500);
         psection = ChaosTools.poincaresos(timeseries, plane; direction=+1, idxs=[2,3]);
         d_traj, v_names = timeseries_to_grid(psection, grid);
         stn, ret_code_stn = create_stn(d_traj, v_names);
         if ret_code_stn == :Success
            break
         end
      end

      P = prob_matrix(stn);
      Q = weight_matrix(stn);

      lyapunov, variance, covariance, ret_code_lyap = analytic_lyapunov(P,Q)
      if ret_code_lyap != :Success
         @show ret_code_lyap, lyapunov, u0
      end

      push!(data0, lyapunov)
      push!(data1, variance)
      push!(data2, covariance)
   
      S, L = network_measures(stn, ensemble, N_steps)
      push!(data3, L)
      push!(data4, S)
end

plot(rho, data3, label = "Random walk")
plot!(rho, data0, label = "Analytic")
plot!(rho, data1, label = "Variance")
plot!(rho, data2, label = "Covariance")
plot!(xlabel = "ρ", ylabel = "Λ", ylim=[-0.5,3.], xlim=[180.,182.], legend = :topleft)


data[data .< 0]
rho[data .< 0]