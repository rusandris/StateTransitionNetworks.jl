using StateTransitionNetworks
using DynamicalSystems
using StatsBase
using Printf

function create_tinkerbell(u0,p)
    tinkerbell_map(u,p,n) = SVector(u[1]^2-u[2]^2+p[1]*u[1]+p[2]*u[2],2*u[1]*u[2]+p[3]*u[1]+p[4]*u[2]) 
    ds = DeterministicIteratedMap(tinkerbell_map,u0,p);
    return ds 
end

function orbit!(orbit_container,ds,nsteps,ntr,u0=initial_state(ds);reinit=false)
    l, = size(orbit_container)
    reinit && reinit!(ds,u0)
    
    l != nsteps && throw(AssertionError("`nsteps`` must be equal to `size(orbit_container)[1]`, Got $(nsteps) and $(l) !"))

	for _ in 1:ntr
		step!(ds)
	end

	for t in 1:nsteps
        step!(ds)

        orbit_container[t,:] .= current_state(ds)[:] 

        #x,y = current_state(ds)
		#orbit_container[t,1] = x
        #orbit_container[t,2] = y
    end

end

#calculate measures for all parameter values in ps 
function calc_measures(ds,ps,orders;out_file_entropy,out_file_lambda,param_index,sts_eltype=Int64,grid_size,grid_edges,u0=initial_state(ds),reinit=false,T,Ttr,ensemble, N_steps,ϵ)

    dim = dimension(ds)
    sts = zeros(sts_eltype,T)
    sts_copy = similar(sts)

    if dim == 1
        orbit = zeros(T)
    else
        orbit = zeros(Float64,T,dim)
    end

    S_vals = zeros(length(orders))
    Λ_vals = similar(S_vals)
    S_sim_vals = similar(S_vals)
    Λ_sim_vals = similar(S_vals)

    #loop through param values 
    for i in 1:length(ps)

        reinit!(ds)
        set_parameter!(ds,param_index,ps[i])

        #step ds to get the orbit
        orbit!(orbit,ds,T,Ttr,u0;reinit=reinit)

        #first order symbolic time series
        timeseries_to_grid!(sts,orbit,grid_size;grid_edges=grid_edges)
        #calculate for every order
       for (j,order) in enumerate(orders)

            if order > 1
                sts_copy .= sts #copy of sts to be mutated by higher_order_symbolics!
                higher_order_symbolics!(sts_copy,order)
                P,_,x = calculate_transition_matrix!(@view sts_copy[1:end-order])
                flush(stderr)
                flush(stdout)
            else
                P,_,x = calculate_transition_matrix!(sts)
            end
            
            if is_stochastic(P) && is_strongly_connected(P)
                S,Λ = network_measures(P;x=x,ϵ=ϵ)
                S_sim,Λ_sim = network_measures(P,ensemble, N_steps)

                S_vals[j] = S
                Λ_vals[j] = Λ
                S_sim_vals[j] = S_sim
                Λ_sim_vals[j] = Λ_sim
            else
                @warn "P is unusable for b = $(ps[i])" * " !"
                flush(stderr)

                S_vals[j] = NaN
                Λ_vals[j] = NaN
                S_sim_vals[j] = NaN
                Λ_sim_vals[j] = NaN
            end

        end

        open(out_file_entropy, "a") do io
            writedlm(io, [ps[i],S_vals...,S_sim_vals...]')
        end  
        
        open(out_file_lambda, "a") do io
            writedlm(io, [ps[i],Λ_vals...,Λ_sim_vals...]')
        end 

        i % 100 == 0 && @info "i = $(i) / $(length(ps)) -> measure calculations for p = $(ps[i]) done."
        flush(stderr)

    end
end

function calc_random_walk_measures(ds,ps,orders;out_file_entropy,out_file_lambda,param_index,sts_eltype=Int64,grid_size,grid_edges,u0=initial_state(ds),reinit=false,T,Ttr,ensemble, N_steps,ϵ)

    dim = dimension(ds)
    sts = zeros(sts_eltype,T)
    sts_copy = similar(sts)

    if dim == 1
        orbit = zeros(T)
    else
        orbit = zeros(Float64,T,dim)
    end

    S_sim_vals = zeros(length(orders))
    Λ_sim_vals = similar(S_sim_vals)

    #loop through param values 
    for i in 1:length(ps)

        reinit!(ds)
        set_parameter!(ds,param_index,ps[i])

        #step ds to get the orbit
        orbit!(orbit,ds,T,Ttr,u0;reinit=reinit)

        #first order symbolic time series
        timeseries_to_grid!(sts,orbit,grid_size;grid_edges=grid_edges)
        #calculate for every order
       for (j,order) in enumerate(orders)

            if order > 1
                sts_copy .= sts #copy of sts to be mutated by higher_order_symbolics!
                higher_order_symbolics!(sts_copy,order)
                P,_,x = calculate_transition_matrix!(@view sts_copy[1:end-order])
                flush(stderr)
                flush(stdout)
            else
                P,_,x = calculate_transition_matrix!(sts)
            end
            
            if is_stochastic(P) && is_strongly_connected(P)
                S_sim,Λ_sim = network_measures(P,ensemble, N_steps)

                S_sim_vals[j] = S_sim
                Λ_sim_vals[j] = Λ_sim
            else
                @warn "P is unusable for b = $(ps[i])" * " !"
                flush(stderr)

                S_sim_vals[j] = NaN
                Λ_sim_vals[j] = NaN
            end

        end

        open(out_file_entropy, "a") do io
            writedlm(io, [ps[i],S_sim_vals...]')
        end  
        
        open(out_file_lambda, "a") do io
            writedlm(io, [ps[i],Λ_sim_vals...]')
        end 

        i % 100 == 0 && @info "i = $(i) / $(length(ps)) -> measure calculations for p = $(ps[i]) done."
        flush(stderr)

    end
end



function calc_lyapunovs(ds,ps,T,Ttr;param_index=1,u0=initial_state(ds))
    #calc lyapunov exponents
    reinit!(ds,u0)
    lyap_exps::Vector{Float64} = zeros(length(ps))

    systems = [deepcopy(ds) for i in 1:Threads.nthreads()]
    pushfirst!(systems, ds) # we can save 1 copy


    Threads.@threads for i in 1:length(ps)
        id = Threads.threadid()
        flush(stderr)
        system = systems[id]
        set_parameter!(system,param_index,ps[i])
        lyap_exps[i] =  lyapunov(system,T;Ttr=Ttr,u0=u0)
    end

    @info "Lyap exponent calc done."
    flush(stderr)
    flush(stdout)
    return lyap_exps
end

function calc_save_renyi_spectrums_grid(ds,qs,ps,orders;param_index,sts_eltype=Int64,grid_size,grid_edges,reinit=false,u0=initial_state(ds),T,Ttr,out_dir="./")

    sts = zeros(sts_eltype,T)
    sts_copy = similar(sts)

    if dimension(ds) == 1
        orbit = zeros(Float64,T)
    else
        orbit = zeros(Float64,T,dimension(ds))
    end
    writedlm(out_dir * "orders" * ".txt",orders)


    #calc spectrums for each parameter value
    for (i,p) in enumerate(ps) 
        reinit!(ds)
        set_parameter!(ds,param_index,p)

        #allocate arrays for spectrums, sts and orbit
        renyi_spectrums = zeros(length(qs),length(orders))
        ps_measures = zeros(2,length(orders))
    
        #step ds to get the orbit
        orbit!(orbit,ds,T,Ttr,u0;reinit=reinit)
        timeseries_to_grid!(sts,orbit,grid_size;grid_edges=grid_edges)   

        #loop through orders
        for j in 1:length(orders)

            if orders[j] > 1
                sts_copy .= deepcopy(sts) #copy of sts to be mutated by higher_order_symbolics!
                higher_order_symbolics!(sts_copy,orders[j])
                P,_,x = calculate_transition_matrix(@view sts_copy[1:end-orders[j]])
            else
                
                P,_,x = calculate_transition_matrix!(sts) #first order P calc remaps sts to smaller symbols
            end

            if is_stochastic(P) && is_strongly_connected(P)
                #writedlm("$(out_dir)" * "symbols_grid_p$(p)" * ".txt",sts)
                renyi_spectrums[:,j] .= renyi_entropy_spectrum(P,qs;x=x,verbose=false)
                ps_measures[:,j] .= network_measures(P;x=x)
            else
                @warn "P not strongly connected or not stochastic for parameter p=$p ..."
                renyi_spectrums[:,j] .= fill(NaN,length(qs))
                ps_measures[:,j] .= fill(NaN,2)
            end

            #@info "Calculation done for order $(orders[j]) ..."
            flush(stderr)
        end
        @info "Spectrum calculations done for parameter p=$p ..."
        flush(stderr)
        T_string = @sprintf "%.E" T
        Ttr_string = @sprintf "%.E" Ttr
        writedlm(out_dir * "renyi_spectrums_grid_param_$(p)" * "_T$T_string" * "_Ttr$Ttr_string" * "_gridsize_$(grid_size)" * ".txt",hcat(qs,renyi_spectrums))
        writedlm(out_dir * "measures_q1_grid_param_$(p)" * "_T$T_string" * "_Ttr$Ttr_string" * "_gridsize_$(grid_size)" * ".txt",ps_measures)
    end
end

function calc_save_renyi_spectrums_op(ds,qs,ps,ws;param_index,sts_eltype=Int64,order=1,u0=initial_state(ds),T,Ttr,out_dir="./")
    dim = dimension(ds)
    orbit = zeros(Float64,T,dim)

    writedlm("$(out_dir)" * "ws_op" * ".txt",ws)

    #calc spectrums for each parameter value
    for (i,p) in enumerate(ps) 
        reinit!(ds)
        set_parameter!(ds,param_index,p)

        #allocate arrays for spectrums, sts and orbit
        renyi_spectrums = zeros(length(qs),length(ws))
        ps_measures = zeros(2,length(ws))
    
        #step ds to get the orbit
        orbit!(orbit,ds,T,Ttr,u0)
           

        #loop through orders
        for (j,w) in enumerate(ws)
       
            sts = codify(OrdinalPatterns{w}(),orbit[:,1])
            sts = Vector{sts_eltype}(sts)

            if order > 1
                higher_order_symbolics!(sts,order) #possible integer symbol overflow 
                P,_,x = calculate_transition_matrix(@view sts[1:end-order+1])
            else
                P,_,x = calculate_transition_matrix!(sts) #first order P calc remaps sts to smaller symbols
            end

            if is_stochastic(P) && is_strongly_connected(P)
            
                #writedlm("$(out_dir)" * "symbols_op_p_$(p)" * "_w_$(w)" * ".txt",sts)

                renyi_spectrums[:,j] .= renyi_entropy_spectrum(P,qs;x=x,verbose=false)
                ps_measures[:,j] .= network_measures(P;x=x)
            else
                @warn "P not strongly connected or not stochastic for parameter p=$p ..."
                renyi_spectrums[:,j] .= fill(NaN,length(qs))
                ps_measures[:,j] .= fill(NaN,2)
            end
            @info "Calculation done for w=$w ..."
            flush(stderr)
        end
        @info "Calculation done for parameter p=$p ..."
        flush(stderr)

        T_string = @sprintf "%.E" T
        Ttr_string = @sprintf "%.E" Ttr
        writedlm(out_dir * "renyi_spectrums_op_param_$(p)" * "_T$T_string" * "_Ttr$Ttr_string" * "_order_$(order)" * ".txt",hcat(qs,renyi_spectrums))
        writedlm(out_dir * "measures_q1_op_param_$(p)" * "_T$T_string" * "_Ttr$Ttr_string" * "_order_$(order)" * ".txt",ps_measures)
    end
end

function walk_length_statistics(P, ensemble::Int64, N_steps::Int64)
    walk_length = zeros(Float64,ensemble)
    for i in 1:ensemble
        walk_length[i] = random_walk_on_stn(P, N_steps)
    end
   	avg = mean(walk_length)/N_steps
    variance = var(walk_length,corrected=false)/N_steps
    return walk_length, avg, variance
end

function walk_length_statistics(ds, ps; param_index=1, sts_eltype=Int64 ,order=1, grid_size,T,u0=initial_state(ds), Ttr, ensemble, N_steps)
    walk_lengths = zeros(ensemble,length(ps))
    means = zeros(length(ps))
    variances = similar(means)
    Ss = similar(means)
    Λs = similar(means)
    sts = zeros(sts_eltype,T+1)

    for (i,p) in enumerate(ps)
        set_parameter!(ds,param_index,p)
        orbit,_ = trajectory(ds,T,u0;Ttr=Ttr)
        sts .= timeseries_to_grid(orbit,grid_size)

        if order > 1
            higher_order_symbolics!(sts,order)
            P,_,x = calculate_transition_matrix(@view sts[1:end-order+1])
        else order == 1
            P,_,x = calculate_transition_matrix(sts)
        end
        
        if is_stochastic(P) && is_strongly_connected(P)
            Ls, avg_norm, variance_norm = walk_length_statistics(P, ensemble, N_steps)
            walk_lengths[:,i] .= Ls
            means[i] = avg_norm
            variances[i] = variance_norm
            S,Λ = network_measures(P;x=x)
            Ss[i] = S
            Λs[i] = Λ
        end
    end
    return walk_lengths,means,variances,Ss,Λs
end
