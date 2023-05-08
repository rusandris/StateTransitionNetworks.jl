include("rotating_plane.jl")

rho_values = [180.1,180.7,180.78,181.1]

parameter_sweep(rho_values;grid_size=20,T=1000,Ttr=500,Δt=0.01,θ_min=0,θ_max=π/2,Δθ=0.001,save = true)



