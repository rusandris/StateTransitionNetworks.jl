using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using GraphPlot

system = Systems.lorenz(ρ=180.70)

# Plane of the Poincare section
plane = (1,15.0);

# Standard Poincare section
psection = poincaresos(
    system, 
    plane, 
    5000; 
    Ttr=500, 
    direction=+1,
    u0=rand(Float64,3).*40 .-20,
    rootkw = (xrtol = 1e-8, atol = 1e-8) 
    );

# Poincare section from samples of the timeseries
traj = trajectory(system, 5000; Δt=0.001, Ttr=500);

psection = poincaresos(
    traj, 
    plane;
    direction=+1,
    );


# Scatter plot of a Poincare section - OPTIONAL
scatter(psection[:, 2], psection[:,3], 
    legend=false, 
    markersize=4, 
    markeropacity=0.1, 
    markerstrokewidth=0,
    color = :blue,
    xlims=(-100,0),
    ylims=(170,210))

# STN creation
grid = 20;
timeseries = psection[:,2:end];
discrete_timeseries, vertex_names = timeseries_to_grid(timeseries, grid);
stn_Q, stn_P = create_STN(discrete_timeseries, vertex_names)

# Obtaining the list of weights
w_list = [];
for e in collect(edges(stn_P))
    push!(w_list,e.weight)
end

# Plotting STN with geometrical data
pl = gplot(
    stn_P,vertex_names[:,2], -vertex_names[:,3],
    nodelabel = 1:nv(stn_P),
    edgelinewidth = w_list .^0.25,
    edgestrokec = "gray",
    NODESIZE = 0.05,
    linetype = "curve",
    nodefillc = "green"
    )