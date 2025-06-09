using DynamicalSystems
using StateTransitionNetworks
#using BenchmarkTools
using Plots,LaTeXStrings

ds = Systems.logistic()
T = Int(1e7)
Ttr= Int(1e4)
grid_edges = [0.0,1.0]
#set_parameter!(ds,1,3.8)
traj,ts = trajectory(ds,T,[0.3];Ttr=Ttr)
x = traj[:,1]

#partition resolutions
ns = [10,20,50,100,300,500,1000,1500,2000,2500,3000]
ms = [2:2:10;]

runtimes_lambda = Float64[]
Ns = Int64[]
nzvals = []
runtimes_sts = Float64[]
runtimes_P = Float64[]

for n in ns 
    @show n
    sts = timeseries_to_grid(traj,n;grid_edges=grid_edges)
    P,_,x = calculate_transition_matrix(sts)
    push!(nzvals,length(P.nzval))
    N = size(P)[1]
    @show N
    push!(Ns,N)
    #_,t,_ = @btimed network_measures(P;x=x)

    #Δt_sts = @belapsed timeseries_to_grid(traj,$n;grid_edges=grid_edges)
    #Δt_P = @belapsed calculate_transition_matrix($sts)
    #Δt_lambda = @belapsed lyapunov_measure($P;x=$x)

    #average over ensemble
    Δt_lambda_avg = 0.0
    Δt_sts_avg = 0.0
    Δt_P_avg = 0.0
    for _ in 1:200
        _,Δt_sts,_ = @timed timeseries_to_grid(traj,n;grid_edges=grid_edges)
        _,Δt_P,_ = @timed calculate_transition_matrix(sts)
        _,Δt_lambda,_ = @timed lyapunov_measure(P;x=x,alg=iterative_linsolve)

        Δt_sts_avg += Δt_sts 
        Δt_P_avg += Δt_P 
        Δt_lambda_avg += Δt_lambda
    end

    #@time lyapunov_measure(P;x=x)
    #push!(runtimes_sts,Δt_sts)
    #push!(runtimes_P,Δt_P)
    push!(runtimes_lambda,Δt_lambda_avg ./ 200)
    push!(runtimes_P,Δt_P_avg ./ 200)
    push!(runtimes_sts,Δt_sts_avg ./ 200)

end

plot!(Ns,nzvals,markershape=:circle,ms=1,yaxis=:log,xaxis=:log,legendfontsize=18,xlabel="N",guidefontsize=20,tickfontsize=15,label=L"r=3.8",ylabel="# trans ")
#plot(Ns,nzvals,markershape=:circle,ms=1,yaxis=:log,xaxis=:log,legendfontsize=18,xlabel="N",guidefontsize=20,tickfontsize=15,label=L"r=4.0",ylabel="# trans ")
plot!(Ns,Ns,markershape=:circle,ms=1,ls=:dash,lc=:black,label=L"N",yaxis=:log,xaxis=:log,xlabel="N",ylabel="")
ylabel!("# trans")
pl = plot!(ylabel="# trans",title="Logistic map",dpi=300)
savefig(pl,"logistic_trans_scaling.png")

plot!(Ns,runtimes_lambda,markershape=:circle,ms=1,yaxis=:log,xaxis=:log,xlabel="N",ylabel="t")
plot!(Ns,runtimes_sts,markershape=:circle,ms=1,yaxis=:log,xaxis=:log)
plot!(Ns,runtimes_P,markershape=:circle,ms=1,yaxis=:log,xaxis=:log)
plot!(Ns,Ns/Ns[1]*runtimes_lambda[1],markershape=:circle,ms=1,ls=:dash,lc=:gray10,yaxis=:log,xaxis=:log)



runtimes_sts = Float64[]
runtimes_P = Float64[]
runtimes_lambda = Float64[]
Ns = Int64[]
nzvals = []
#order
sts = timeseries_to_grid(traj,5;grid_edges=grid_edges)
for m in ms 
    @show m
    sts_copy = deepcopy(sts)
    higher_order_symbolics!(sts_copy,m)
    P,_,x = calculate_transition_matrix(@view sts_copy[1:end-m+1])
    push!(nzvals,length(P.nzval))
    N = size(P)[1]
    @show N
    push!(Ns,N)
    #_,t,_ = @btimed network_measures(P;x=x)

    #Δt_sts = @belapsed timeseries_to_grid(traj,$n;grid_edges=grid_edges)
    #Δt_P = @belapsed calculate_transition_matrix($sts)
    #Δt_lambda = @belapsed lyapunov_measure($P;x=$x)

    #average over ensemble
    #average over ensemble
    Δt_lambda_avg = 0.0
    Δt_sts_avg = 0.0
    Δt_P_avg = 0.0
    for _ in 1:50
        _,Δt_sts,_ = @timed timeseries_to_grid(traj,5;grid_edges=grid_edges)
        _,Δt_P,_ = @timed calculate_transition_matrix(sts)
        _,Δt_lambda,_ = @timed lyapunov_measure(P;x=x,alg=iterative_linsolve)

        Δt_sts_avg += Δt_sts 
        Δt_P_avg += Δt_P 
        Δt_lambda_avg += Δt_lambda
    end

    push!(runtimes_lambda,Δt_lambda_avg ./ 50)
    push!(runtimes_P,Δt_P_avg ./ 50)
    push!(runtimes_sts,Δt_sts_avg ./ 50)
end

plot!(Ns,runtimes_lambda,markershape=:circle,lw=2,ms=1,xlabel=L"N(m)",label=L"\Lambda"*" calc",guidefontsize=20,legendfontsize=15,ylabel="runtime")
plot!(Ns,runtimes_sts,markershape=:circle,ms=1,lw=2,xlabel=L"N(m)",label="binning",guidefontsize=20,legendfontsize=15)
plot!(Ns,runtimes_P,markershape=:circle,ms=1,lw=2,xlabel=L"N(m)",label="trans matrix",guidefontsize=20,legendfontsize=15)

plot!(dpi=300)
plot!(Ns,Ns/Ns[1]*runtimes_lambda[1],label="",markershape=:circle,ms=1,ls=:dash,lc=:gray10)

savefig(pl,"logistic_lyapunov_runtime_scaling.png")
plot!(ms,runtimes_lambda,markershape=:circle,ms=1,yaxis=:log,xlabel="m",ylabel="t")
plot!(ms,10 ,markershape=:circle,ms=1,ls=:dash,lc=:gray10,yaxis=:log)
