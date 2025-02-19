using SparseArrays
function early_warning_signals(timeseries,w,lags;τ=1)
    ts = timeseries.parent[timeseries.indices[1]]
    stn_measures = complexity_measures(timeseries,w,τ)

    #autocorrelation and variance
    v = var(ts;corrected=false)
    #lag1, normalized by variance
    #substracts mean
    ac = autocor(timeseries,[1])[1]
    acfs = acf_sum(timeseries,lags)

    return [stn_measures...,v,ac,acfs]
end

function early_warning_signals(timeseries,grid_size,grid_edges,lags;transform=false,outside_grid=:error,nbins=50,force_strongly_connected=false)

    ts = timeseries.parent[timeseries.indices[1]]
    if transform
        ts = transform_timeseries(ts;nbins=nbins)
    end
    stn_measures = complexity_measures(ts,grid_size,grid_edges;force_strongly_connected=force_strongly_connected,outside_grid=outside_grid)

    #autocorrelation and variance
    v = var(ts;corrected=false)
    #lag1, normalized by variance
    #substracts mean
    ac = autocor(ts,[1])[1]
    acfs = acf_sum(ts,lags)

    return [stn_measures...,v,ac,acfs]
end

acf_sum(x,lags) = sum(autocor(x,lags))

function transform_timeseries(timeseries;nbins=50,GCDF=false,standardize=false)

    #standardize
    if standardize==true
        standardize!(timeseries)
    end

    if GCDF == false
        #create histogram
        h = fit(Histogram, vec(timeseries), nbins=nbins) 
        h_normed = normalize(h;mode=:probability)

        #get bin properties, probabilities
        bin_width = Float64(h_normed.edges[1].step)
        bin_edges = h_normed.edges[1]
        bin_probs = h_normed.weights
        
        #calculate cumulative distribution
        cumulative_probs = cumsum(bin_probs)
        #pushfirst!(cumulative_probs,0.0)
        push!(cumulative_probs,1.0)

        #interpolator and transform to uniformly distributed timeseries
        interp = LinearInterpolator(bin_edges,cumulative_probs)
        timeseries_transformed = interp.(timeseries)
    else
        timeseries_transformed = normcdf.(0.0,1.0,timeseries)
    end

    return timeseries_transformed,cumulative_probs,bin_edges 
end

#return measures of a timeseries using grid
function complexity_measures(timeseries,grid_size::Int,grid_edges::Vector{Float64};force_strongly_connected=true,size_ratio_threshold=0.9,outside_grid=:error)

	symbolic_timeseries = timeseries_to_grid(timeseries,grid_size;grid_edges=grid_edges,outside_grid=outside_grid)
    unique_symbols = unique(symbolic_timeseries)
    nr_unique_symbols = length(unique_symbols) 
	P,_,ρ_est = calculate_transition_matrix(symbolic_timeseries)
    nr_trans = length(P.nzval)
	#Λ = lyapunov_measure(P;x=ρ_est)[1]

    S,Λ,C1,C2 = complexity_measures(P,ρ_est;size_ratio_threshold = size_ratio_threshold,force_strongly_connected=force_strongly_connected)
    return S,Λ,C1,C2,nr_unique_symbols,nr_trans
end

#return measures from a symbolic timeseries
function complexity_measures(symbolic_timeseries;force_strongly_connected=true,size_ratio_threshold = 0.9)
    unique_symbols = unique(symbolic_timeseries)
    nr_unique_symbols = length(unique_symbols) 
	P,_,ρ_est = calculate_transition_matrix(symbolic_timeseries)
    nr_trans = length(P.nzval)
	#Λ = lyapunov_measure(P;x=ρ_est)[1]
    S,Λ,C1,C2 = complexity_measures(P::SparseMatrixCSC,ρ_est;force_strongly_connected=force_strongly_connected,size_ratio_threshold = size_ratio_threshold)
    return S,Λ,C1,C2,nr_unique_symbols,nr_trans
end

#return measures from a timeseries using OP
function complexity_measures(timeseries,w::Int,τ::Int;force_strongly_connected=false,size_ratio_threshold=0.9)

    symbolic_timeseries = codify(OrdinalPatterns{w}(τ),timeseries[:,1])
    unique_symbols = unique(symbolic_timeseries)
    nr_unique_symbols = length(unique_symbols) 
    P,_,ρ_est = calculate_transition_matrix(symbolic_timeseries)
    nr_trans = length(P.nzval)

    S,Λ,C1,C2 = complexity_measures(P,ρ_est;size_ratio_threshold = size_ratio_threshold,force_strongly_connected=force_strongly_connected)

    return S,Λ,C1,C2,nr_unique_symbols,nr_trans
end

#return measures from a timeseries using transition matrix P (internal) 
function complexity_measures(P::SparseMatrixCSC,ρ_est;force_strongly_connected=true,size_ratio_threshold = 0.9)
    n0 = size(P)[1]

    #checking stochasticity and connectedness
    #if not stochastic => not connected
     
    if is_stochastic(P) && is_strongly_connected(P)
        S,Λ =  network_measures(P;x=ρ_est)
        C1,C2 = bit_number_measures(ρ_est)

    elseif force_strongly_connected
        P = make_strongly_connected(P)
        n = size(P)[1]

        #if strongly connected component is too small, return NaNs instead
        if n/n0 < size_ratio_threshold
            println("size ratio is $(n/n0)" * " < $size_ratio_threshold")
            S,Λ,C1,C2 = (NaN,NaN,NaN,NaN)
        else
            ρ = stationary_distribution(P)
            S,Λ =  network_measures(P;x=ρ)
            C1,C2 = bit_number_measures(ρ)
        end

    else
        S,Λ,C1,C2 = (NaN,NaN,NaN,NaN)
    end
    return S,Λ,C1,C2
end

#partition resolution sensitivity analysis

function prsa(timeseries,part_resolution_range = Vector{Int64};other_part_param)
    measures = zeros(length(part_resolution_range),4)

    for (i,res) in enumerate(part_resolution_range)
        measures[i,:] .= complexity_measures(timeseries,res,other_part_param)
    end

    return measures
end

function standardize!(data::Matrix)
    dims = size(data)[2]
    l = size(data)[1]
    means,stds = mean_and_std(data,1) 

    for d in 1:dims
        for i in 1:l
            data[i,d] = (data[i,d] - means[d])/stds[d]
        end
    end

end

function standardize!(data)
    l = length(data)
    mean,std = mean_and_std(data,1) 
    data .= (data .- mean) ./ std
end

function translate_mean!(data::Matrix)
    dims = size(data)[2]
    l = size(data)[1]
    means = mean(data;dims=1) 

    for d in 1:dims
        for i in 1:l
            data[i,d] = (data[i,d] - means[d])
        end
    end

end