using StateTransitionNetworks
using Plots
using DynamicalSystems
using SparseArrays
using StatsBase

using LinearAlgebra
using DelimitedFiles
using LaTeXStrings
using Printf
cd(@__DIR__)

figs_dir="../figs/"
mkpath(figs_dir)

data_dir="../data/supplimentary/crit_tent_maps/"
mkpath(data_dir)

include("functions_chaotic_maps.jl")

#-----------------------------------critical map--------------------------------------------

function f_crit(u, p, n)
    r = p[1]
    x = u[1]
    return SVector(1 - abs(x^r-(1-x)^r)^(1.0/r))
end


function static_renyi_entropy(ρ, q)
    Hq = 0.
    if isapprox(q,1.0)
        for i in 1:length(ρ)
            Hq += -ρ[i]*log(ρ[i])
        end
        return Hq
    else
        for i in 1:length(ρ)
            Hq += (ρ[i])^q
        end
        Hq = log(Hq)/(1-q)
        return Hq
    end
end

function static_renyi_entropy_spectrum(ρ, qs=collect(0:0.01:2))
    Hqs = zeros(length(qs))
    for (i,q) in enumerate(qs)
        Hqs[i] = static_renyi_entropy(ρ, q)
    end
    return Hqs
end

###############################
### Measures for a single value
###############################

r = 2.
ds_crit = DeterministicIteratedMap(f_crit, [0.4], [r])


const T::Int64 = 10^7
const T_string::String = @sprintf "%.E" T
const Ttr::Int64 = 10^5
const Ttr_string::String = @sprintf "%.E" Ttr
const grid_size::Int64 = 2^5
const grid_edges::Vector{Float64} = [0,1]
writedlm(data_dir*"method_params.txt",[T,Ttr,grid_size])

timeseries, = trajectory(ds_crit, T, [0.3]; Ttr=Ttr);
#=
# plot time series
plot(timeseries[:,1][end-500:end], label=L"r=%$(r)", linewidth=2, alpha=0.8, color=:red)
plot!(xlabel=L"\tau", ylabel=L"x_\tau", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, legendfontsize=16, dpi=300)
plot!(ylim=[0,1])
savefig("$(path)/critical_x-t_r=$(r)_tmax=$(T_string)_ttrans=$(Ttr_string).pdf")
savefig("$(path)/critical_x-t_r=$(r)_tmax=$(T_string)_ttrans=$(Ttr_string).svg")


sts = timeseries_to_grid(timeseries, grid_size; grid_edges = [0., 1.]);

# 1st order
order = 1
higher_order_symbolics!(sts, order)
@time P,Q,x = calculate_transition_matrix(@view sts[1:end-order+1]);
@time S, Λ = network_measures(P; x=x, ϵ=1e-2, maxiter=10^8, alg=hybrid_solve)
@time C1, C2 = bit_number_measures(x) 

# histogram
histogram(timeseries[:,1]; bins=grid_size, normalize=:pdf, color=:red, alpha=0.8, lw=2, la=1, label="")
plot!(xlabel=L"x", ylabel=L"\rho(x)", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, legendfontsize=16, dpi=300)
savefig("$(path)/critical_natural_dist_r=$(r)_tmax=$(T_string)_ttrans=$(Ttr_string).svg")

#@time P,Q,x = calculate_transition_matrix_no_remap(sts);
#ε = 1/grid_size
#plot([ε/2+(ε)*i for i in 1:grid_size], x/ε)
#plot!(ylim=[0,1])

# higher order
order = 10
sts2 = zeros(Int128, length(sts))
sts2 .= sts
@time higher_order_symbolics!(sts2, order)
nr_symbols = length(unique(sts2))
@time P,Q,x = calculate_transition_matrix(@view sts2[1:end-order+1]);
@time S10, Λ10 = network_measures(P; x=x, ϵ=1e-2, maxiter=10^8, alg=hybrid_solve)

λ = lyapunov(ds_crit, T; Ttr=Ttr, d0 = 1e-9)

f_name = path*"/critical_renyi-q_dq=0.01_r=$(r)_tmax=10^8_ttrans=10^6_grid=$(grid_size)_order=10.dat"
data = readdlm(f_name)
qs = data[:,1]
Hs10 = data[:,2]
#i1 = findall(qs .≈ 1.)
Hs10[1]

println("Critical map")
println("state entropy: ", C1)
println("state varentropy: ", C2)
println("Lyapunov exponent: ", λ)
println("KS entropy (1st order): ", S)
println("KS entropy (10th order): ", S10)
println("Topological entropy (10th order): ", Hs10[1])
println("Lyapunov measure (1st order): ", Λ)
println("Lyapunov measure (10th order: ", Λ10)
=#

# compute DYNAMICAL Renyi entropy spectrum for HIGHER ORDER states
grid_size = 2^5
qs = collect(0.01:0.01:2)
Os = [1,2,4,8,10]
sts = zeros(Int128,length(timeseries))
timeseries_to_grid!(sts, timeseries, grid_size; grid_edges = [0., 1.]);
sts2 = zeros(Int128, length(sts))
for (i,o) in enumerate(Os[1:end])
        sts2 .= sts
        @show o
        @time higher_order_symbolics!(sts2, o)
        nr_symbols = length(unique(sts2))
        @show nr_symbols
        @time P,Q,x = calculate_transition_matrix(@view sts2[1:end-o+1]);
        Hs = renyi_entropy_spectrum(P, qs; x=x, verbose=true)
        f_name = data_dir*"critical_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(Int(o)).dat"
        writedlm(f_name,hcat(qs, Hs))
end


#######################################
# compute STATIC Renyi entropy spectrum
#######################################
grid_size = 2^5
ε = 1/grid_size
qs = collect(0.01:0.1:2)
#Hqs = [ComplexityMeasures.entropy(Renyi(; q=q, base=ℯ), ValueHistogram(ε), timeseries) for q ∈ qs]
plot(qs,Hqs, label=L"n=%$(grid_size)", linewidth=3, color="red")
plot!(xlabel=L"q", ylabel=L"H_q", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, legendfontsize=16, dpi=300)
plot!(xlim=[0,2], ylim=[0.,7.], legend=:bottomleft, left_margin=3Plots.mm)
#annotate!(pl, (-0.18,  1), text("(a)", :left, 24))
savefig(figs_dir*"critical_static_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size).pdf")


sts = timeseries_to_grid(timeseries, grid_size; grid_edges = [0., 1.]);
@time P,Q,x = calculate_transition_matrix(sts);
Hqs = static_renyi_entropy_spectrum(x, qs)
plot!(qs,Hqs)

# higher order states
grid_size = 2^5
qs = collect(0.01:0.01:2)
Os = [1,2,4,8,10]
sts = zeros(Int128,length(timeseries))
timeseries_to_grid!(sts, timeseries, grid_size; grid_edges = [0., 1.]);
sts2 = zeros(Int128, length(sts))
for (i,o) in enumerate(Os[1:end])
        sts2 .= sts
        @show o
        @time higher_order_symbolics!(sts2, o)
        nr_symbols = length(unique(sts2))
        @show nr_symbols
        @time P,Q,x = calculate_transition_matrix(@view sts2[1:end-o+1]);
        Hs = static_renyi_entropy_spectrum(x, qs)
        f_name = data_dir*"critical_static_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(Int(o)).dat"
        writedlm(f_name,hcat(qs, Hs))
end



#=
###################
# fractal dimension
###################
timeseries, = trajectory(ds_crit, 10^7, [0.3]; Ttr=10^6);

ες = 2 .^ (-15:0.5:5) # semi-random guess
generalized_dim(timeseries, ες; q=1)

function generalized_entropy_spectrum(X, qs)
    Ds = zeros(length(qs))
    for (i,q) in enumerate(qs)
        @show q
        Ds[i] = generalized_dim(X, ες; q=q)
    end
    return Ds

end

qs = collect(0.01:0.1:2)
Ds = generalized_entropy_spectrum(timeseries, qs)
plot!(qs, Ds)


# correlation sum
timeseries, = trajectory(ds_crit, 10^4, [0.3]; Ttr=10^6);
ες = 2 .^ (-15:0.5:5) # semi-random guess

Cs = correlationsum(timeseries, ες; show_progress=true)

xs = log2.(ες)
ys = log2.(Cs)
scatter(xs, ys; ylabel = L"\log(C_2)", xlabel = L"\log (\epsilon)")
Δ = slopefit(xs, ys, LargestLinearRegion())
grassberger_proccacia_dim(timeseries, estimate_boxsizes(timeseries); )
grassberger_proccacia_dim(timeseries, ες; )


function generalized_correlation_sum_spectrum(X, qs)
    Cs = zeros(length(qs))
    for (i,q) in enumerate(qs)
        @show q
        #Cs[i] = correlationsum(X, ες; q=q)
        Cs[i] = grassberger_proccacia_dim(X, ες; q=q)

    end
    return Cs
end
qs = collect(0.01:0.1:2)
Cs = generalized_correlation_sum_spectrum(timeseries, qs)
plot!(qs,Cs)
=#

#-----------------------------------tent map--------------------------------------------

function f_tent(u, p, n)
    r = p[1]
    x = u[1]
    if x<=r
        return SVector(x/r)
    else
        return SVector((1-x)/(1-r))
    end
end

function analytical_renyi_entropy(w, q)
    if q==0
        return log(2)
    elseif q==1
        return -w*log(w)-(1-w)*log(1-w)
    else
        return log(w^q+(1-w)^q)/(1-q)
    end
end

function analytical_renyi_entropy_spectrum(w, qs; verbose=true)
    Hs = zeros(length(qs))
    for (i,q) in enumerate(qs)
        verbose && @show q
        Hs[i] = analytical_renyi_entropy(w, q)
    end
    return Hs
end

r = 0.8
ds_tent = DeterministicIteratedMap(f_tent, [0.4], [r])

timeseries, = trajectory(ds_tent, T, [0.3]; Ttr=Ttr);
#=
# time series plot
plot(timeseries[:,1][end-500:end], label=L"r=%$(r)", linewidth=2, alpha=0.8, color=:orange)
plot!(xlabel=L"\tau", ylabel=L"x_\tau", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, legendfontsize=16, dpi=300)
savefig(figs_dir*"asymmetric_triangular_x-t_r=$(r)_tmax=$(T_string)_ttrans=$(Ttr_string).pdf")
savefig(figs_dir*"asymmetric_triangular_x-t_r=$(r)_tmax=$(T_string)_ttrans=$(Ttr_string).svg")

# 1st order
order = 1
higher_order_symbolics!(sts, order)
@time P,Q,x = calculate_transition_matrix(@view sts[1:end-order+1]);
@time S, Λ = network_measures(P; x=x, ϵ=1e-2, maxiter=10^8, alg=hybrid_solve)
@time C1, C2 = bit_number_measures(x) 
=#

sts = timeseries_to_grid(timeseries, grid_size; grid_edges = [0., 1.]);

# compute Renyi entropy spectrum for HIGHER ORDER states
grid_size = 2^5
qs = collect(0.01:0.01:2)
Os = [1,2,4,8,10]
writedlm(data_dir*"orders_tent_dynamic.txt",Os)
sts = zeros(Int128,length(timeseries))
timeseries_to_grid!(sts, timeseries, grid_size; grid_edges = [0., 1.]);
sts2 = zeros(Int128, length(sts))
for (i,o) in enumerate(Os[1:end])
        sts2 .= sts
        @show o
        @time higher_order_symbolics!(sts2, o)
        nr_symbols = length(unique(sts2))
        @show nr_symbols
        @time P,Q,x = calculate_transition_matrix(@view sts2[1:end-o+1]);
        Hs = renyi_entropy_spectrum(P, qs; x=x, verbose=true)
        f_name = data_dir*"asymmetric_triangular_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(Int(o)).dat"
        writedlm(f_name,hcat(qs, Hs))
end


#######################################
# compute STATIC Renyi entropy spectrum
#######################################
grid_size = 2^5
ε = 1/grid_size
qs = collect(0.01:0.1:2)
#=
#Hqs = [ComplexityMeasures.entropy(Renyi(; q=q, base=ℯ), ValueHistogram(ε), timeseries) for q ∈ qs]
plot(qs,Hqs, label=L"n=%$(grid_size)", linewidth=3, color="red")
plot!(xlabel=L"q", ylabel=L"H_q", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, legendfontsize=16, dpi=300)
plot!(xlim=[0,2], ylim=[0.,7.], legend=:bottomleft, left_margin=3Plots.mm)
#annotate!(pl, (-0.18,  1), text("(a)", :left, 24))
savefig(figs_dir*"asymmetric_triangular_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size).pdf")
=#

sts = timeseries_to_grid(timeseries, grid_size; grid_edges = [0., 1.]);

# compute for higher order states
grid_size = 2^5
qs = collect(0.01:0.01:2)
Os = [8,10,12,14] #[1,2,4,8,10]
writedlm(data_dir*"orders_tent_static.txt",Os)
sts = zeros(Int128,length(timeseries))
timeseries_to_grid!(sts, timeseries, grid_size; grid_edges = [0., 1.]);
sts2 = zeros(Int128, length(sts))
for (i,o) in enumerate(Os[1:end])
        sts2 .= sts
        @show o
        @time higher_order_symbolics!(sts2, o)
        nr_symbols = length(unique(sts2))
        @show nr_symbols
        @time P,Q,x = calculate_transition_matrix(@view sts2[1:end-o+1]);
        Hs = static_renyi_entropy_spectrum(x, qs)
        f_name = data_dir*"asymmetric_triangular_static_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(Int(o)).dat"
        writedlm(f_name,hcat(qs, Hs))
end

