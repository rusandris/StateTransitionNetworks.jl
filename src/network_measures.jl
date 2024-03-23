export random_walk_on_stn, randomwalk_step,network_measures, sinai_kolmogorov_entropy, measure_convergence, lyapunov_measure, stationary_distribution
export bit_number_measures,renyi_entropy,renyi_entropy_spectrum


"""
	network_measures(P::AbstractMatrix) -> S, Λ
Calculates the Sinai-Kolmogorov Entropy (S) and Lyapunov measure (Λ) from P by using the analytical definitions of both quantities.

Optional Keyword arguments:
* `x`: probability distribution of states. If nothing is provided, it is calculated from `P'*x = x` using `stationary_distribution` (`Krylovkit`'s `eigsolve`).
* `alg`: linear solver (`hybrid_solve` by default) for Λ -> `iterative_linsolve`, `linsolve` (`KrylovKit`),`hybrid_solve` are the options
* `ϵ`: tolerance for Λ calculation
* `maxiter`: maxiter for Λ calculation
"""
function network_measures(P::AbstractMatrix;x=nothing,ϵ=1e-12,maxiter=1000,alg=hybrid_solve)
   	entropy, ret_code_entr = sinai_kolmogorov_entropy(P;x=x)
    lyapunov, variance, covariance, ret_code_lyap = lyapunov_measure(P;x=x,ϵ=ϵ,maxiter=maxiter,alg=alg)
	return entropy, lyapunov
end

"""
	sinai_kolmogorov_entropy(stn::MetaDiGraph) -> S
Calculates analytically the Sinai-Kolmogorov entropy on a STN by extracting the Q weight matrix and P transition probability matrix.
"""
function sinai_kolmogorov_entropy(stn)
    Q = get_weight_matrix(stn)
    P = get_transition_matrix(stn)

    entropy = -sum(Q .* -sparse_log(P))
    return entropy
end

"""
	stationary_distribution(P) -> x
Calculates the stationary probability based on the probability matrix P. Each element of the resulting vector is the probability of finding the system in node i.
"""
function stationary_distribution(P::AbstractMatrix)
	vals, vecs, info = eigsolve(P',1,:LR)
	info.converged < 1 && @warn "KrylovKit.eigsolve did not converge!" 
	x = real.(vecs[1]) ./sum(real.(vecs[1]))
	return x
end

"""
	sinai_kolmogorov_entropy(Q,P) -> S
Calculates analytically the Sinai-Kolmogorov entropy given the Q weight matrix and P transition probability matrix of the STN. 

"""
function sinai_kolmogorov_entropy(Q,P)
    entropy = -sum(Q .* -sparse_log(P))
	return entropy
end

"""
	sinai_kolmogorov_entropy(P) -> S
Calculates analytically the Sinai-Kolmogorov entropy given the P transition probability matrix of the STN. 

"""
function sinai_kolmogorov_entropy(P::AbstractMatrix;x=nothing)
		
	if isnothing(x)
		x = stationary_distribution(P)
	end
	
	v = ones(length(x))

	L = -sparse_log(P)
	L = P.*L
	entropy = (x'*L*v)[1]
	return real(entropy), :Success
end


"""
	lyapunov_measure(P) -> Λ
Calculates analytically the Lyapunov measure given the the P transition probability matrix. Uses `iterative_linsolve to find needed matrix inverses.`

Optional Keyword arguments:
* `x`: probability distribution of states. If nothing is provided, it is calculated from `P'*x = x` using `stationary_distribution` (`Krylovkit`'s `eigsolve`).
* `alg`: linear solver (`hybrid_solve` by default) -> `iterative_linsolve`, `linsolve` (`KrylovKit`),`hybrid_solve` are the options
* `ϵ`: tolerance for linear solvers
* `maxiter`: maxiter for linear solvers
"""
function lyapunov_measure(P::AbstractMatrix;x=nothing,alg=hybrid_solve,ϵ=1e-12,maxiter=1000)
	
	if isnothing(x)
		x = stationary_distribution(P)
	end
	
	xt = transpose(x)
	v = ones(length(x))
	
	L = -sparse_log(P)
	L2 = P.*L.^2
	L = P.*L
	
	X = PseudoDenseMatrix(x) 
	
	if alg == iterative_linsolve || alg == hybrid_solve
		z,convergence_info = alg(P,X,L*v;ϵ = ϵ,maxiter=maxiter)
		
		if convergence_info isa KrylovKit.ConvergenceInfo
			convergence_info.converged != 1 && @warn "KrylovKit.linsolve did not converge!"
		elseif convergence_info isa Bool
			convergence_info || @warn "iterative_linsolve did not converge! Your system converges slowly/might not converge at all. Try setting `maxiter` kwarg to bigger or `ϵ` to a higher value!"
		end
		
	elseif alg == KrylovKit.linsolve
		z,info = KrylovKit.linsolve(I - P + X,L*v)
		info.converged != 1 && @warn "KrylovKit.linsolve did not converge!" 
	else
		error("This algorithm is not implemented yet! See the function's docstring for available options.")
	end
	
	covariance = 2*xt*L*(z - X*L*v)
	
	variance = (xt*L2*v)[1] - (xt*L*v)[1]^2
	lyapunov =  variance + covariance[1]
	if imag(covariance) < 1.0e-3
	   return lyapunov, variance, covariance, :Success
	else
	   @show covariance
	   return real(lyapunov), real(variance), real(covariance), :ComplexCovariancveWarning
	end
end

"""
	sparse_log(A::AbstractMatrix) -> log(A)
Calculates the piecewise logarithm of a sparse matrix without explicitely converting it to dense matrix.
"""
function sparse_log(A::AbstractMatrix)
	P = copy(A)
	vals = nonzeros(P)
	m, n = size(P)
	for j = 1:n
		for i in nzrange(P, j)
			vals[i] = log(vals[i])
		end
	end
	return P
end
 
 
#Schlogl notation: C1,C2

function bit_number_measures(x::Vector{Float64})

	l = -log.(x)
	replace!(l, Inf=>0.0)
	entropy = sum((x .* l)) #C1
    variance =  sum(x .* l .* l) - entropy^2 #C2
	
	return real(entropy), real(variance), :Success
end

#----------------Renyi entropy spectrum---------------

function renyi_entropy(P::SparseMatrixCSC{Float64, Int64}, q::Float64; x::Vector{Float64}=stationary_distribution(P),tol::Float64=1e-8, maxiter::Int64=10^4,verbosity::Int64=0)

	P_q = P.^q
	l, = eigsolve(P_q; verbosity=verbosity, issymmetric=false, ishermitian=false, tol=tol, maxiter=maxiter)

	λ_max, = findmax(abs.(l))

	# for PPC there are multiple eigenvalues of the same magnitude and some are complex
	#abs(imag(l[i_max]))<0.01 ? H = log(λ_max)/(1-q) : H = -1

	H = log(abs(λ_max))/(1-q)
	return H
	
end


function renyi_entropy(P::SparseMatrixCSC{Float64, Int64}, q::Float64,n::Float64; x::Vector{Float64}=stationary_distribution(P))

	x_q = x .^ q
	P_q = P .^ q
	v_q = ones(length(x_q))
	return log((x_q' * P_q^n * v_q))/(n*(1-q))

end

function renyi_entropy_spectrum(P::SparseMatrixCSC{Float64, Int64}, qs::Vector{Float64}; x::Vector{Float64}=stationary_distribution(P), verbose=true)
    Hs = zeros(length(qs))
    for (i,q) in enumerate(qs)
        verbose && @show q
		if q == 1.0
			Hs[i] = sinai_kolmogorov_entropy(P;x=x)[1]
			continue 
		end
        Hs[i] = renyi_entropy(P, q; x=x)
    end
    return Hs
end
