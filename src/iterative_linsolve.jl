export iterative_linsolve, hybrid_solve, PseudoDenseMatrix



#--------------pseudo dense matrices-------------------
#IDEA: to avoid storing all the elements of the vx
#outer product, store only a vector 

struct PseudoDenseMatrix
	v::Vector{Float64}
end

function Base.:*(M::PseudoDenseMatrix,S::SparseMatrixCSC)
	return PseudoDenseMatrix(transpose(transpose(M.v)*S))
end

function Base.:*(M::PseudoDenseMatrix,v::Vector{Float64})
	return transpose(M.v)*v * ones(length(v))
end

function Base.:*(v::Vector{Float64},M::PseudoDenseMatrix)
	return transpose(v)*M.v
end

function Base.:+(S::SparseMatrixCSC,M::PseudoDenseMatrix)
	return Matrix(S .+ transpose(M.v))
end

Base.:+(M::PseudoDenseMatrix,S::SparseMatrixCSC) = S + M


"""
	iterative_linsolve(S::SparseMatrixCSC,D::PseudoDenseMatrix,y::AbstractVector;z0::AbstractVector = rand(length(y)), ϵ = 1e-12,maxiter=1000) -> z::AbstractVector
Solves for a linear system of equations of the form `(S+D)z = y` iteratively. Here `S` is a sparse matrix, `D` is a dense matrix with rank 1 (PseudoDenseMatrix).
Not guaranteed to converge.
"""
function iterative_linsolve(S::SparseMatrixCSC,D::PseudoDenseMatrix,y::Vector{Float64},z0::Vector{Float64} = rand(length(y)); ϵ = 1e-12,maxiter=1000)
	
	converged = false
	norm_diff = Inf
	z = z0
	i = 0
	
	while i < maxiter
		z = y + S*z - D*z

		norm_diff = norm(z-z0)
		norm_diff/length(z) > 1.0  && error("Norm of z-z0/length(z) is bigger than 1!")
		#@show norm_diff/length(z)
		
		i+=1
		
		if norm_diff < ϵ 
			converged = true
			return z,converged
		end
		z0 = z
	end
	return z,converged
end

"""
	hybrid_solve(S::SparseMatrixCSC,D::PseudoDenseMatrix,y::AbstractVector;z0::AbstractVector = rand(length(y)), ϵ = 1e-12,maxiter=1000) -> z::AbstractVector
Solves for a linear system of equations of the form `(S+D)z = y`. Here `S` is a sparse matrix, `D` is a dense matrix with rank 1 (PseudoDenseMatrix). If size of S is bigger than `Ns x Ns` than it uses `iterative_linsolve`, otherwise `KrylovKit.linsolve` is used
"""
function hybrid_solve(S::SparseMatrixCSC,D::PseudoDenseMatrix,y::AbstractVector,z0::AbstractVector = rand(length(y)); ϵ = 1e-12,maxiter=1000,Ns=1000,fallback=true)
	N = size(S)[1]
	
	#if size of matrix N is small, try KrylovKit.linsolve
	if N < Ns
		return linsolve(I - S + D,y;tol = ϵ,maxiter=maxiter)
	
	#else try own iterative_linsolve
	else
		z, converged = iterative_linsolve(S,D,y;ϵ = ϵ,maxiter=maxiter)
		
		#if that didn't converge go back to fallback case (KrylovKit.linsolve)
		if !converged && fallback
			return linsolve(I - S + D,y;tol = ϵ,maxiter=maxiter)
		else
			return z, converged
		end
	end
end
