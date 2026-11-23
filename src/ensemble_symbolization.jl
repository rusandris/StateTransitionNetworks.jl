export ensemble_measures

# for parameter range using BINS with :grid or OP with :OP (returns vectors of measures for each parameter value)
function ensemble_measures(parameter_range::Vector{Float64},N_ens::Int,create_system,symbolization_pars::NamedTuple,order::Int;
	init_ensemble_grid_edges,Ttr=100,init_sys_pars=[],p_index=1,Δt=1)
    
    # unwrap symbolization_pars
    par_names = keys(symbolization_pars)
    symb_type = symbolization_pars.symb_type
    symbolization_settings = symbolization_pars[filter(p-> p !=:symb_type ,par_names)]
    
    @info "Calculating ensemble measures for parameter interval = [$(parameter_range[1]), $(parameter_range[end])] with $(N_ens) ensemble members and order = $(order) using $(symb_type) symbolization"
    
    Ss = Float64[]
    Λs = Float64[]
    C1s = Float64[]
    C2s = Float64[]
    # initialize the ensemble system
    u0 = zeros(N_ens, 2);
    ds_ens = Array{Any}(undef, N_ens);

    init_ics!(u0, N_ens, init_ensemble_grid_edges)
    create_ensemble!(ds_ens, u0, create_system; N_ens=N_ens, p_index=p_index, p_value=parameter_range[1], init_sys_pars=init_sys_pars)

    # transient steps until the attractor is reached
    @info "Computing transient steps until the attractor is reached for Ttr = $(Ttr) with Δt = $(Δt)..."
    step_ensemble!(ds_ens, Int(floor(Ttr/Δt)); Δt=Δt, verbose=false)
    
    @showprogress for i in 1:length(parameter_range)
        @info "Setting parameter index $(p_index) to value $(parameter_range[i]) for all ensemble members..."
        # set the new parameter value
        set_parameter!.(ds_ens, p_index, parameter_range[i])
        # more transients due to the parameter change
        step_ensemble!(ds_ens, Int(floor(Ttr/Δt)); Δt=Δt)
        
        # compute measures for a single step from HIGHER order symbols
        symb_ens = zeros(Int128, N_ens, order)
        ho_symb_ens = zeros(Int128, N_ens, 2)
        
        if symb_type == :grid
            symbolize_ensemble!(ds_ens, symb_ens, ho_symb_ens, order; N_ens=length(ds_ens), Δt=Δt, symbolization_settings...)
            symb_stat = ensemble_transition_statistics!(ds_ens, symb_ens, ho_symb_ens, order;
            	N_ens=length(symb_ens[:,1]), Δt=Δt, symbolization_settings...)
        elseif symb_type == :OP
        	w = symbolization_settings.w
        	τ = symbolization_settings.τ
            ts_ens = timeseries_ens(ds_ens, ((w-1)*τ+order)*Δt; Δt=Δt, exact_times=true, Ttr=0.0, var_idx=1)
            # first 1:(w-1)*τ+order times for initialization
            symbolize_ensemble!(ts_ens[:,1:(w-1)*τ+order], symb_ens, ho_symb_ens, order; w=w, τ=τ, Δt=Δt)
            # second 2:(w-1)*τ+order+1 times for next symbolic step
            symb_stat = ensemble_transition_statistics!(ts_ens[:,order+1:(w-1)*τ+order+1], symb_ens, ho_symb_ens, order; Δt=Δt,symbolization_settings...)
        else
        	throw(ArgumentError("Invalid symbolization type, expected :grid or :OP, got $(symb_type)"))    
        end
              
        P = calculate_transition_matrix(symb_stat)
        ρ = calculate_symbol_probabilities(symb_stat)
        if !is_strongly_connected(P)
            @info "Transition matrix is not strongly connected, making it strongly connected..."
            P = make_strongly_connected(P)
            ρ = stationary_distribution(P)
        end
        S, Λ = network_measures(P; x=ρ, ϵ=1e-2, maxiter=10^8, alg=hybrid_solve)
        C1, C2 = bit_number_measures(ρ) 
        push!(Ss,S)
        push!(Λs,Λ)
        push!(C1s,C1)
        push!(C2s,C2)
    end
    return Ss, Λs, C1s, C2s
end

# method for nonauto case using :grid or :OP with (returns vectors of measures for each step)
function ensemble_measures(n_steps::Int,N_ens::Int,create_system,symbolization_pars::NamedTuple,order::Int;
	init_ensemble_grid_edges,Ttr=100,auto_sys_pars=[],nonauto_sys_pars,p_index=1,Δt=1)
    
    # unwrap symbolization_pars
    par_names = keys(symbolization_pars)
    symb_type = symbolization_pars.symb_type
    symbolization_settings = symbolization_pars[filter(p-> p !=:symb_type,par_names)]
    
    Ss = Float64[]
    Λs = Float64[]
    C1s = Float64[]
    C2s = Float64[]
    # initialize the ensemble system
    u0 = zeros(N_ens, 2);
    ds_ens = Array{Any}(undef, N_ens);

    init_ics!(u0, N_ens, init_ensemble_grid_edges)
    create_ensemble!(ds_ens, u0, create_system; N_ens=N_ens, init_sys_pars=auto_sys_pars)

    # transient steps until the attractor is reached
    @info "Computing transient steps until the attractor is reached for Ttr = $(Ttr) with Δt = $(Δt)..."
    step_ensemble!(ds_ens, Int(floor(Ttr/Δt)); Δt=Δt, verbose=false)
    set_time!.(ds_ens, 0.0)
    
    states = current_state.(ds_ens)
    @show extrema(states[1])
    
    #set to nonauto (drift)
    for i in 1:length(ds_ens)
    	set_parameters!(ds_ens[i], nonauto_sys_pars)
    end 
        
	# init symbolized ensemble with HIGHER order symbols
    symb_ens = zeros(Int128, N_ens, order)
    ho_symb_ens = zeros(Int128, N_ens, 2)
    
    if symb_type == :grid
            symbolize_ensemble!(ds_ens, symb_ens, ho_symb_ens, order; N_ens=length(ds_ens), Δt=Δt, symbolization_settings...)
    elseif symb_type == :OP
        	w = symbolization_settings.w
        	τ = symbolization_settings.τ
            ts_ens = timeseries_ens(ds_ens, ((w-1)*τ+order)*Δt; Δt=Δt, exact_times=true, Ttr=0.0, var_idx=1)
            # first 1:(w-1)*τ+order times for initialization
            symbolize_ensemble!(ts_ens[:,1:(w-1)*τ+order], symb_ens, ho_symb_ens, order; w=w, τ=τ, Δt=Δt)     
    end
    
    @showprogress for n in 0:n_steps
        #ε = ds_ens[1].integ.p[3] + ds_ens[1].integ.p[4]*n*2π
    	#@show n, ε

        if symb_type == :grid
            symb_stat = ensemble_transition_statistics!(ds_ens, symb_ens, ho_symb_ens, order;
            	N_ens=length(symb_ens[:,1]), Δt=Δt, symbolization_settings...)
        elseif symb_type == :OP
        	w = symbolization_settings.w
        	τ = symbolization_settings.τ
            ts_ens = timeseries_ens(ds_ens, ((w-1)*τ+order)*Δt; Δt=Δt, exact_times=true, Ttr=0.0, var_idx=1)
            symb_stat = ensemble_transition_statistics!(ts_ens[:,order+1:(w-1)*τ+order+1], symb_ens, ho_symb_ens, order; Δt=Δt,symbolization_settings...)
        else
        	throw(ArgumentError("Invalid symbolization type, expected :grid or :OP, got $(symb_type)"))    
        end
              
        P = calculate_transition_matrix(symb_stat)
        ρ = calculate_symbol_probabilities(symb_stat)
        if !is_strongly_connected(P)
            @info "Transition matrix is not strongly connected, making it strongly connected..."
            P = make_strongly_connected(P)
            ρ = stationary_distribution(P)
        end
        S, Λ = network_measures(P; x=ρ, ϵ=1e-2, maxiter=10^8, alg=hybrid_solve)
        C1, C2 = bit_number_measures(ρ) 
        push!(Ss,S)
        push!(Λs,Λ)
        push!(C1s,C1)
        push!(C2s,C2)
    end
    return Ss, Λs, C1s, C2s
end

#=
# NEW FUNCTION for a single parameter value using BINS with :grid or OP with :OP options (returns a single value for each measure)
function ensemble_measures(parameter_value::Float64,N_ens::Int,create_system,symbolization_pars::Tuple,order::Int;
	init_ensemble_grid_edges,Ttr=100,init_sys_pars=[],p_index=1,Δt=1)
    # unwrap symbolization_pars
    symb_type = symbolization_pars[1]
    if symb_type == :grid
        grid_size = symbolization_pars[2]
        grid_edges = symbolization_pars[3]
    elseif symb_type == :OP
        w = symbolization_pars[2]
        τ = symbolization_pars[3]
    else
        throw(ArgumentError("Invalid symbolization type, expected :grid or :OP, got $(symb_type)"))    
    end
    @info "Calculating ensemble measures for parameter value = $(parameter_value) with $(N_ens) ensemble members and order = $(order) using $(symb_type) symbolization"

    # initialize the ensemble system
    u0 = zeros(N_ens, 2);
    ds_ens = Array{Any}(undef, N_ens);

    init_ics!(u0, N_ens, init_ensemble_grid_edges)
    init_ensemble!(ds_ens, u0, create_system; N_ens=N_ens, p_index=p_index, p_value=parameter_value, init_sys_pars=init_sys_pars)

    # transient steps until the attractor is reached
    @info "Computing transient steps until the attractor is reached for Ttr = $(Ttr) with Δt = $(Δt)..."
    step_ensemble!(ds_ens, Int(floor(Ttr/Δt)); Δt=Δt)
    
    # compute measures for a single step from HIGHER order symbols
    symb_ens = zeros(Int128, N_ens, order)
    ho_symb_ens = zeros(Int128, N_ens, 2)

    if symb_type == :grid
        symbolize_ensemble!(ds_ens, symb_ens, ho_symb_ens, order; N_ens=length(ds_ens), grid_size=grid_size, grid_edges=grid_edges, Δt=Δt)
        symb_stat = ensemble_transition_statistics!(ds_ens, symb_ens, ho_symb_ens, order; N_ens=length(symb_ens[:,1]), grid_size=grid_size, grid_edges=grid_edges, Δt=Δt)
    elseif symb_type == :OP
        ts_ens = timeseries_ens(ds_ens, ((w-1)*τ+order)*Δt; Δt=Δt, exact_times=true, Ttr=0.0, var_idx=1)
        # first 1:(w-1)*τ+order times for initialization
        symbolize_ensemble!(ts_ens[:,1:(w-1)*τ+order], symb_ens, ho_symb_ens, order; w=w, τ=τ, Δt=Δt)
        # second 2:(w-1)*τ+order+1 times for next symbolic step
        symb_stat = ensemble_transition_statistics!(ts_ens[:,order+1:(w-1)*τ+order+1], symb_ens, ho_symb_ens, order; w=w, τ=τ, Δt=Δt)
    end 

    P = calculate_transition_matrix(symb_stat)
    ρ = calculate_symbol_probabilities(symb_stat)
    if !is_strongly_connected(P)
        @info "Transition matrix is not strongly connected, making it strongly connected..."
        P = make_strongly_connected(P)
        ρ = stationary_distribution(P)
    end
    S, Λ = network_measures(P; x=ρ, ϵ=1e-2, maxiter=10^8, alg=hybrid_solve)
    C1, C2 = bit_number_measures(ρ) 
    return S, Λ, C1, C2
end
=#
# generate a trajectory for each ensemble member and return a matrix of time series (N_ens x length of the time series) for a given variable index
function timeseries_ens(ds_ens, T; Δt = 0.1, exact_times = true, Ttr = 0.0, var_idx=1)
    N_ens = length(ds_ens)
    t0 = current_time(ds_ens[1])
    tvec = (t0 + Ttr):Δt:(t0 + T + Ttr)
    ts_ens = Matrix{Float64}(undef, N_ens, length(tvec))
    if Ttr > 0
        step_ens!(ds_ens, Int(floor(Ttr/Δt)); Δt=Δt, exact_times=exact_times)
    end
    for (i, t) in enumerate(tvec)
        step_ensemble!(ds_ens; Δt=Δt, exact_times=exact_times)
        for j in 1:N_ens
            ts_ens[j,i] = current_state(ds_ens[j])[var_idx]
        end
        if !all(successful_step.(ds_ens))
            # Diverged trajectory; set final state to remaining set
            for k in (i + 1):length(tvec)
                ts_ens[:,k] .= ts_ens[:,i]
            end
            break
        end
    end

    return ts_ens
end

# map to symbolvia paritioning the state space into a grid of size `grid_size` and mapping the state to a symbol based on which grid cell it falls into
function map_to_symbol(u; grid_size=20, grid_edges=[-1.5,-1.5,1.5,1.5], method=error_outside_grid)
    if length(grid_edges) != 4
		throw(ArgumentError("Length of `grid_edges` should be 4 in case of 2D time series, got $(length(grid_edges)) !"))
	end
    dim = Int(length(grid_edges)/2)
	grid_steps = zeros(dim)
    for i in 1:dim 
        max = grid_edges[i+dim] #pick out max of the grid along a given axis
        min = grid_edges[i] #pick out min of the grid along a given axis

        #add padding to right border of the grid
        #so the points on the edge aren't singled out 
        #maxs_plus[i] = nextfloat(max,grid_border_iter*nfloat)
        #grid cell width
        grid_steps[i] = (max - min)/grid_size
        #test if largest symbol goes out of bounds
        #max_symbols[i] = assign_grid_cell(max,min,grid_steps[i])
    end
    x_min,y_min,x_max,y_max = grid_edges
    j_cell = method(u[1],[x_min,x_max],grid_size,grid_steps[1])
    i_cell = method(u[2],[y_min,y_max],grid_size,grid_steps[2])
    symbol = (j_cell-1)*grid_size + i_cell
    return symbol 
end

# map to OP symbols for a given window size w and time delay τ
function map_to_symbol(us::Vector{Float64}; w=3, τ=1)
    if length(us) != ((w-1)*τ+1)
        throw(ArgumentError("Length of `us` should be equal to `(w-1)*τ+1`, got $(length(us)) and $((w-1)*τ+1)!"))
    end
    symbol = codify(OrdinalPatterns{w}(τ),us)[1]
    return symbol
end

# initialize the initial conditions for the ensemble of systems on a grid in the state space defined by `grid_edges` (a vector of tuples of min and max values for each dimension)
function init_ics!(u0, N_ens, grid_edges::Vector{Tuple{Float64, Float64}})
    x0_min,x0_max = grid_edges[1]
    v0_min,v0_max = grid_edges[2]
    N_ens_1D = Int(sqrt(N_ens))
    for i in 1:N_ens_1D
        v0 = v0_min + (v0_max-v0_min)/N_ens_1D*(i-1)
        for j in 1:N_ens_1D  
            x0 = x0_min + (x0_max-x0_min)/N_ens_1D*(j-1)
            k = (i-1)*N_ens_1D+j
            u0[k,:] .= x0, v0
        end
    end    
end

# initialize the ensemble of systems with initial conditions `u0` and set the parameter at index `p_index` to `p_value` for each system in the ensemble
function create_ensemble!(ds_ens, u0, create_system; N_ens=length(ds_ens), p_index=0, p_value=0.0, init_sys_pars=[])
    for i in 1:N_ens
        if init_sys_pars != []
            ds_ens[i] = create_system(u0[i,:]; init_sys_pars=init_sys_pars)
        else
            ds_ens[i] = create_system(u0[i,:])
        end
        p_index != 0 && set_parameter!(ds_ens[i], p_index, p_value)
        set_time!(ds_ens[i], 0.0)
    end
end

# set the time of the ensemble of systems to `t`
function set_time!(ds, t=0.0)
    ds.integ.t = t
end

# step the ensemble of systems for `n` steps with time step `Δt` and optionally print the current step if `verbose` is true
function step_ensemble!(ds_ens, n=1; Δt=1, exact_times=true, verbose=false)
    for step in 1:n
        if n>1 && verbose
            @show step
        end
        step!.(ds_ens, Δt, exact_times)
    end
end

# compute transition statistics for the ensemble using higher order symbols for a given order, grid size and grid edges
function ensemble_transition_statistics!(ds_ens, symb_ens, ho_symb_ens, order; N_ens=length(ds_ens), grid_size=20, grid_edges=[-1.5,-1.5,1.5,1.5], Δt=1, return_probabilities=false)
    symbol_space_size = grid_size^2
    ho_symbol_space_size = get_ho_symb_space_size(symbol_space_size, order)
    f_prev = zeros(ho_symbol_space_size)
    f = zeros(ho_symbol_space_size)
    symb_stat = SymbolStatistics(typeof(symb_ens[1]), typemax(Int64))
    for i in 1:N_ens
        # move the HIGHER order symbols in the buffer to the left
        ho_symb_ens[i,1] = ho_symb_ens[i,2]
        # move the FIRST order symbols in their buffer to the left
        symb_ens[i,1:order-1] .= symb_ens[i,2:order]
        # step the integrator and compute the next symbol
        step!(ds_ens[i], Δt, true)
        s = map_to_symbol(current_state(ds_ens[i]); grid_size=grid_size, grid_edges=grid_edges)
        # add new FIRST order symbol to the buffer
        symb_ens[i,order] = s
        # compute new higher order symbol
        ho_symb_ens[i,2] = higher_order_state(symb_ens[i,1:order], symbol_space_size; warn_overflow=true)[1]
        # add transition to the SymbolStatistics container
        ho_sprev = ho_symb_ens[i,1]
        ho_s = ho_symb_ens[i,2]
        add_transition!(symb_stat, (ho_sprev,ho_s))
        f_prev[ho_sprev] += 1
        f[ho_s] += 1
    end
    if return_probabilities
        return symb_stat, f_prev, f
    else
        return symb_stat
    end
end

# compute the maximum number of higher order symbols for a given symbol space size and order
function get_ho_symb_space_size(symb_space_size, order)
    ho_symb_space_size = 0
    for i in 1:order
        ho_symb_space_size += symb_space_size^i
    end
    return ho_symb_space_size
end


# compute transition statistics for the ensemble using higher order symbols for a given order, window size w and time delay τ
function ensemble_transition_statistics!(ts_ens::Matrix{Float64}, symb_ens, ho_symb_ens, order; w, τ, Δt=1, return_probabilities=false)
    # a SINGLE symbol is computed from a window of the time series, so the length of the time series should be (w-1)*τ+1
    if length(ts_ens[1,:]) != ((w-1)*τ+1)
        throw(ArgumentError("Length of `ts` should be equal to `(w-1)*τ+1`, got $(length(ts_ens[1,:])) and $((w-1)*τ+1)!"))
    end
    N_ens=size(ts_ens,1)
    symbol_space_size = factorial(w)
    ho_symbol_space_size = get_ho_symb_space_size(symbol_space_size, order)
    f_prev = zeros(ho_symbol_space_size)
    f = zeros(ho_symbol_space_size)
    symb_stat = SymbolStatistics(typeof(symb_ens[1]), typemax(Int64))
    for i in 1:N_ens
        # move the HIGHER order symbols in the buffer to the left
        ho_symb_ens[i,1] = ho_symb_ens[i,2]
        # move the FIRST order symbols in their buffer to the left
        symb_ens[i,1:order-1] .= symb_ens[i,2:order]
        # symbolize the next segment in time seriesn
        s = map_to_symbol(ts_ens[i,1:(w-1)*τ+1], w=w, τ=τ)
        # add new FIRST order symbol to the buffer
        symb_ens[i,order] = s
        # compute new higher order symbol
        ho_symb_ens[i,2] = higher_order_state(symb_ens[i,1:order], symbol_space_size; warn_overflow=true)[1]
        # add transition to the SymbolStatistics container
        ho_sprev = ho_symb_ens[i,1]
        ho_s = ho_symb_ens[i,2]
        add_transition!(symb_stat, (ho_sprev,ho_s))
        f_prev[ho_sprev] += 1
        f[ho_s] += 1
    end
    if return_probabilities
        return symb_stat, f_prev, f
    else
        return symb_stat
    end
end

# symbolize the ensemble using higher order symbols for a given order, grid size and grid edges
function symbolize_ensemble!(ds_ens, symb_ens, ho_symb_ens, order; N_ens=length(ds_ens), grid_size=20, grid_edges=[-1.5,-1.5,1.5,1.5], Δt=1)
    symbol_space_size = grid_size^2
    symb_ens[:,1] .= map_to_symbol.(current_state.(ds_ens); grid_size=grid_size, grid_edges=grid_edges)
    for j in 2:order
        step_ensemble!(ds_ens; Δt=Δt)
        symb_ens[:,j] .= map_to_symbol.(current_state.(ds_ens); grid_size=grid_size, grid_edges=grid_edges)
    end
    for i in 1:N_ens
        ho_symb_ens[i,2] = higher_order_state(symb_ens[i,1:order], symbol_space_size; warn_overflow=true)[1]
    end
end

# compute transition statistics for the ensemble using higher order symbols for a given order, window size w and time delay τ
function symbolize_ensemble!(ts_ens::Matrix{Float64}, symb_ens, ho_symb_ens, order; w, τ, Δt=1)
    # ORDER number of symbols are computed from a window of the time series, so the length of the time series should be (w-1)*τ+order
     if length(ts_ens[1,:]) != ((w-1)*τ+order)
        throw(ArgumentError("Length of `ts` should be equal to `(w-1)*τ+order`, got $(length(ts_ens[1,:])) and $((w-1)*τ+order)!"))
    end
    N_ens=size(ts_ens,1)
    symbol_space_size = factorial(w)
    #ho_symbol_space_size = (symbol_space_size)^order
    
    for j in 1:order
        for i in 1:N_ens
            # symbolize the next segment in time series
            s = map_to_symbol(ts_ens[i,j:(w-1)*τ+j], w=w, τ=τ)
            # add new FIRST order symbol to the buffer
            symb_ens[i,j] = s
        end
    end
    for i in 1:N_ens
        ho_symb_ens[i,2] = higher_order_state(symb_ens[i,1:order], symbol_space_size; warn_overflow=true)[1]
    end
end


