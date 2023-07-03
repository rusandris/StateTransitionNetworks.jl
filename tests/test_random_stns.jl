using StateTransitionNetworks
using Plots
using DynamicalSystems
using Graphs
using GraphPlot
using Random
using Statistics

using LinearAlgebra
using DelimitedFiles
using LaTeXStrings

N = 4
ens = 1000


function random_transition_matrix(N, A)
    P = zeros(N,N)
    for i in 1:N
        for j in 1:N
            if A[i,j]==1
                P[i,j] = rand()
            end
        end
    end
    for i in 1:N
        if sum(P[i,:]) !=0 
            P[i,:] = P[i,:]/sum(P[i,:])
        end
    end
    return P
end

a = 2^(N*N)-1
A = Float64.(reshape(digits(a, base=2, pad=N*N), (N,N)))
P = random_transition_matrix(N,A)
network_measures(P)
stn, retcode = create_stn(P; make_ergodic=true)
P = prob_matrix(stn)
network_measures(P)

data0 = [];
data1 = [];
data2 = [];
data3 = [];
data4 = [];
for a in 1:2^(N*N)-1
    @show a
    A = Float64.(reshape(digits(a, base=2, pad=N*N), (N,N)))
    P = random_transition_matrix(N,A)
    stn, retcode = create_stn(P; make_ergodic=true)
    if retcode == :Success
        av_entropy = 0.
        av_lyapunov = 0.
        var_entropy = 0.
        var_lyapunov = 0.
        for i in 1:ens
            P = random_transition_matrix(N,A)
            stn, retcode = create_stn(P; make_ergodic=true)
            P = prob_matrix(stn)
            entropy, lyapunov = network_measures(P)
            av_entropy += entropy
            av_lyapunov += lyapunov
            var_entropy += entropy^2
            var_lyapunov += lyapunov^2
        end
    else
        av_entropy = -1. * ens
        av_lyapunov = -1. * ens
        var_entropy = 1. * ens
        var_lyapunov = 1. * ens
    end
    push!(data0, a)
    push!(data1, av_entropy/ens)
    push!(data2, av_lyapunov/ens)
    push!(data3, sqrt(var_entropy/ens-(av_entropy/ens)^2))
    push!(data4, sqrt(var_lyapunov/ens-(av_lyapunov/ens)^2))
end

mask1 = data1 .!= -1
plot(data0[mask1], data3[mask1], label=L"\mathrm{std}(S)")
plot!(data0[mask1], data1[mask1], label=L"\langle S \rangle")
plot!(xlabel=L"a", ylabel=L"S", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, xtickfontsize=12, lw=2, legendfontsize=16)
savefig("random_networks_entropy_N=4_uniform.pdf")
mask2 = data2 .!= -1
plot(data0[mask2], data4[mask2], label=L"\mathrm{std}(\Lambda)")
plot!(data0[mask2], data2[mask2], label=L"\langle \Lambda \rangle")
plot!(ylim=[0.01,2])
plot!(xlabel=L"a", ylabel=L"\Lambda", xguidefontsize=22, yguidefontsize=22, tickfontsize=14, xtickfontsize=12, lw=2, legendfontsize=16)
savefig("random_networks_lyapunov_N=4_uniform.pdf")


idx = findall(x->x>1.15, data1)
A = reshape(digits(idx[1], base=2, pad=N*N), (N,N))
P = random_transition_matrix(N,A)
stn, retcode = create_stn(P; make_ergodic=true)
P = prob_matrix(stn)
network_measures(P)
plot_stn(stn)

idx = findall(x->x>10, data4)

# 0 lyapunov measure
idx = findall(x->x≈0, data2)
A = reshape(digits(idx[end], base=2, pad=N*N), (N,N))
P = random_transition_matrix(N,A)
stn, retcode = create_stn(P; make_ergodic=true)
P = prob_matrix(stn)
network_measures(P)

idx = findmax(data1)[end]
A = reshape(digits(idx[end], base=2, pad=N*N), (N,N))
P = random_transition_matrix(N,A)
stn, retcode = create_stn(P; make_ergodic=true)
P = prob_matrix(stn)
network_measures(P)