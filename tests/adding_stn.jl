using Revise
using StateTransitionNetworks
using DynamicalSystems
using ChaosTools
using Plots
using LinearAlgebra
using StatsBase
using MetaGraphsNext
using Graphs

import StateTransitionNetworks.lyapunov_measure
import SparseArrays: spzeros 
import StateTransitionNetworks.calculate_transition_matrix
import Graphs: DiGraph


function timeseries_to_common_grid(timeseries, grid, x_min, x_max, y_min, y_max)    
    T = length(timeseries[:,1])
    dx = 0.5*(x_max-x_min)/grid
    dy = 0.5*(y_max-y_min)/grid
    x_grid = range(x_min-dx, x_max+dx, grid);
    x_min = x_grid[1]
    y_grid = range(y_min-dy, y_max+dy, grid);
    y_min = y_grid[1]
    x_n = Vector{Int64}(undef, T)
    y_n = Vector{Int64}(undef, T)
    x_n, y_n = [], [];

    for row in eachrow(timeseries)
        y = floor(Int,(row[2]-y_min)/Float64(y_grid.step))+1
        x = floor(Int,(row[1]-x_min)/Float64(x_grid.step))+1
        push!(x_n, x)
        push!(y_n, y)
    end
    d_timeseries = [x_n y_n]
    return d_timeseries
end



function add_timeseries(dt_list, grid;make_ergodic=false,verbose=false)
    Q_null = zeros(Int32, grid*grid, grid*grid) # Null matrix with all possible transitions
    vertex_names = [];  # Future name of vertices
    vertex_place = [];  # Rows, and column number in Q_null for a given vertex
    M = zeros(grid,grid)
    nr_vertices = 0;
    for dt in dt_list
        states = []
        O = zeros(Int32, grid*grid, grid*grid)
        for row in eachrow(dt)
            v = (row[1]-1)*grid + row[2]
            x, y = row
            if M[y,x] == 0
                nr_vertices += 1 
                push!(vertex_names, [nr_vertices, x, y])
                push!(vertex_place, v)
                M[y,x] = 1
            end    
            push!(states, v) 
        end
        #count transitions
        next_states = circshift(states,-1)
        for i in eachindex(states[1:end-1])
            O[states[i],next_states[i]] += 1
        end
        Q_null = Q_null + O
    end
    vertex_names = reduce(hcat, vertex_names)'
    Q = Q_null[vertex_place, vertex_place]
	#weight and transition probability matrices
    P = spzeros(nr_vertices, nr_vertices)
	#normalize Q and fill P by normalizing rows
    Q = Q./sum(Q)

	P = calculate_transition_matrix(Q)
	#create directed metagraph with static label and metadata types and default weight 0
	stn, ret_code = create_stn(P;make_ergodic=false,verbose=false)
    for v in 1:nr_vertices
		x,y = vertex_names[v,2:end] 
		stn[v] = Dict(:x => x,:y => y)
	end
    return stn, ret_code
end

function get_grid_edges(psections)
    x_max, y_max = maximum(reduce(vcat, maximum.(Matrix.(psections); dims=1)); dims=1)
    x_min, y_min = minimum(reduce(vcat, minimum.(Matrix.(psections); dims=1)); dims=1)
    return x_min, x_max, y_min, y_max
end

Δt = 0.001;
plane = (1,15.0);
grid = 20;
T = 500;
ρ = 180.1;

u0 = rand(Float64,3).*50 .-25';
system = PredefinedDynamicalSystems.lorenz(u0; ρ=ρ);
timeseries_1 = trajectory(system, T; Δt=Δt, Ttr=1000);
psection_1 = DynamicalSystemsBase.poincaresos(timeseries_1, plane; direction=+1, save_idxs=[2,3]);
dtraj_1, v1 = timeseries_to_grid(psection_1, grid);
stn1, ret_code_stn = create_stn(dtraj_1, v1; make_ergodic=true,verbose=false);

u0 = rand(Float64,3).*50 .-25';
system = PredefinedDynamicalSystems.lorenz(u0; ρ=ρ);
timeseries_2 = trajectory(system, T; Δt=Δt, Ttr=1000);
psection_2 = DynamicalSystemsBase.poincaresos(timeseries_2, plane; direction=+1, save_idxs=[2,3]);
dtraj_2, v2 = timeseries_to_grid(psection_2, grid);
stn2, ret_code_stn = create_stn(dtraj_2, v2; make_ergodic=true,verbose=false);


psections = [psection_1, psection_2];
x_min, x_max, y_min, y_max = get_grid_edges(psections)

dt_1 = timeseries_to_common_grid(psection_1, grid, x_min, x_max, y_min, y_max);
dt_2 = timeseries_to_common_grid(psection_2, grid, x_min, x_max, y_min, y_max);
stn_added1, ret_code = add_timeseries([dt_1, dt_2], grid; make_ergodic=true, verbose=false);
stn_added2, ret_code = add_timeseries([dt_2, dt_1], grid; make_ergodic=true, verbose=false);

plot_stn(stn1;filename="stn1.pdf",nodesize=1,nodefillc="orange",linetype="curve",max_edgelinewidth=1, nodelabels=true)
plot_stn(stn2;filename="stn2.pdf",nodesize=1,nodefillc="green",linetype="curve",max_edgelinewidth=1, nodelabels=true)
plot_stn(stn_added1;filename="stn_added.pdf",nodesize=1,nodefillc="red",linetype="curve",max_edgelinewidth=1, nodelabels=true)
plot_stn(stn_added2;filename="stn_added2.pdf",nodesize=1,nodefillc="salmon",linetype="curve",max_edgelinewidth=1, nodelabels=true)

Q = weight_matrix(stn1)
Q = weight_matrix(stn2)
Q1 = weight_matrix(stn_added1)
Q2 = weight_matrix(stn_added2)
P1 = prob_matrix(stn_added1)
lyapunov_measure(P1)
P2 = prob_matrix(stn_added2)
lyapunov_measure(P2)

Q1.nzval
Q1[2, :]
Q1[2, :].nzval
[Q1[i,:].nzval for i in 1:nv(stn1)]
reduce(vcat, [P1[i,:].nzval for i in 1:nv(stn1)])
minimum(Q2.nzval)^0.125