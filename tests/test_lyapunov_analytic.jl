using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using GraphPlot
using Random
using Statistics

using LinearAlgebra
using DelimitedFiles

function uniquetol(A; kws...)
   S = []
   for a in A
        if !any(s -> isapprox(s, a; kws...), S)
            push!(S, a)
        end
   end
   return S
end

# Lyapunov measure for a given ρ value
Δt = 0.001;
plane = (1,15.0);
grid = 20;
ρ = 180.9;
u0 = rand(3)
T = 5000;
system = Systems.lorenz(u0; ρ=ρ);
timeseries,  = trajectory(system, T; Δt=Δt, Ttr=500);

psection = DynamicalSystemsBase.poincaresos(timeseries, plane; direction=+1, save_idxs=[2,3]);
d_traj, v_names = timeseries_to_grid(psection, grid);
stn, ret_code = create_stn(d_traj, v_names);
ret_code
#plot_stn(stn)

P = prob_matrix(stn);
Q = weight_matrix(stn);

# OLD function for the analytic Lyapunov===============================================================

function analytic_lyapunov_old(P,Q)
   λ, v = eigen(Matrix(P))
   
   if length(uniquetol(λ))<size(P)[1]
      @show λ
   end
   
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
   
   #=
   for (i,λi) in enumerate(λ)
      if abs(λi) ≈ 1
         #@show λi
         Ω[i,i] = 0.
      end
   end
   =#
   
   L = Matrix(-log.(P))
   replace!(L, Inf=>0.0)
   L2 = P.*L.^2
   L = P.*L

   covariance = x*L*v*Ω*inv(v)*L*ones(length(λ))
   #avg = sinai_kolmogorov_entropy(Q,P)
   #sq_avg = sum(Q[Q .!=0] .* (log.(P[P .!=0])).^2)
   #variance = sq_avg - avg^2
   variance = x*L2*ones(length(λ))-(x*L*ones(length(λ)))^2
   lyapunov =  variance + 2*real(covariance)
   if imag(covariance) < 1.0e-3
      return lyapunov, variance, real(covariance), :Success
   else
      @show covariance
      return lyapunov, variance, real(covariance), :ComplexCovariancveWarning
   end
end

# Analytic calculation of the Lyapunov measure ====================================================

function analytic_lyapunov(P)
   λ, V = eigen(Matrix(P))
   λ, X = eigen(transpose(Matrix(P)))
   
   if real(λ[end]) ≈ 1
      x = transpose(X[:,end]./sum(X[:,end]))
      v = V[:,end]./V[1,end]
   else
      return -1, -1, -1, :StochasticMatrixError
   end

   L = Matrix(-log.(P))
   replace!(L, Inf=>0.0)
   L2 = P.*L.^2
   L = P.*L
   I = Diagonal(ones(length(x)))
   S = (I-v*x)+(P-v*x)*inv(I-P+v*x)
   covariance = x*L*S*L*v
   variance = x*L2*v-(x*L*v)^2
   lyapunov =  variance + 2*real(covariance)
   if imag(covariance) < 1.0e-3
      return lyapunov, variance, real(covariance), :Success
   else
      @show covariance
      return lyapunov, variance, real(covariance), :ComplexCovariancveWarning
   end
end


analytic_lyapunov_old(P,Q)
analytic_lyapunov(P)
network_measures(stn, 1000, 10^4)[2]


# Analytic Lyapunov measure over an interval for the Lorenz system ======================================================================================

rho = 180.:0.001:182.;
T = 100000;
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
         timeseries,  = trajectory(system, T; Δt=Δt, Ttr=500);
         psection = DynamicalSystemsBase.poincaresos(timeseries, plane; direction=+1, save_idxs=[2,3]);
         d_traj, v_names = timeseries_to_grid(psection, grid);
         stn, ret_code_stn = create_stn(d_traj, v_names);
         if ret_code_stn == :Success
            break
         end
      end

      P = prob_matrix(stn);
      Q = weight_matrix(stn);

      lyapunov, variance, covariance, ret_code_lyap = analytic_lyapunov(P)
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
plot!(rho, 2 .* data2, label = "2 x Covariance")
plot!(xlabel = "ρ", ylabel = "Λ", ylim=[-1.,3.], xlim=[180.,182.], legend = :topleft)
#plot!(xlabel = "ρ", ylabel = "Λ", ylim=[-0.5,3.], xlim=[180.85,181.], legend = :topleft)

data = hcat(rho, data3, data0, data1, data2)
writedlm("lorenz_analytic_lyapunov_rho=180-182_drho=0.001_t=10^5.dat",data)

# Analytic Lyapunov measure over an interval for the Henon system ======================================================================================

a_values = 1:0.001:1.4;
b = 0.3;
T = 10000;
ensemble = 10000;
N_steps = 10000;
data0 = [];
data1 = [];
data2 = [];
data3 = [];
data4 = [];
for a in a_values
      @show a
      while true # STN must be 'healthy'
         #u0 = rand(Float64,2);
         system = Systems.henon([0.,0.]; a=a, b=b);
         timeseries,  = trajectory(system, T; Ttr=1000);
         d_traj, v_names = timeseries_to_grid(timeseries, grid);
         stn, ret_code_stn = create_stn(d_traj, v_names);
         if ret_code_stn == :Success
            break
         end
      end

      P = prob_matrix(stn);
      Q = weight_matrix(stn);

      lyapunov, variance, covariance, ret_code_lyap = analytic_lyapunov(P)
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

plot(a_values, data3, label = "Random walk")
plot!(a_values, data0, label = "Analytic")
plot!(a_values, data1, label = "Variance")
plot!(a_values, 2 .* data2, label = "2 x Covariance")
plot!(xlabel = "a", ylabel = "Λ", ylim=[-0.25,1.], xlim=[1.,1.4], legend = :topleft)

plot(a_values, data4, label = "Entropy")
plot!(a_values, data1, label = "Variance")

data = hcat(a_values, data3, data0, data1, data2, data4)
writedlm("henon_analytic_lyapunov_a=1.0-1.4_da=0.001_t=10^5.dat",data)