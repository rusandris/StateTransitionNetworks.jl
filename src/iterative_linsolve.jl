export iterative_linsolve, PseudoDenseMatrix


#--------------pseudo dense matrices-------------------
#IDEA: to avoid storing all the elements of the vx
#outer product, store only a vector 

struct PseudoDenseMatrix
	v::AbstractVector
end

function Base.:*(M::PseudoDenseMatrix,S::SparseMatrixCSC)
	return PseudoDenseMatrix(transpose(transpose(M.v)*S))
end

function Base.:*(M::PseudoDenseMatrix,v::AbstractVector)
	return transpose(M.v)*v * ones(length(v))
end

function Base.:*(v::AbstractVector,M::PseudoDenseMatrix)
	return transpose(v)*M.v
end

function Base.:+(S::SparseMatrixCSC,M::PseudoDenseMatrix)
	return Matrix(S .+ transpose(M.v))
end

Base.:+(M::PseudoDenseMatrix,S::SparseMatrixCSC) = S + M


"""
	iterative_linsolve(S::SparseMatrixCSC,D::PseudoDenseMatrix,y::AbstractVector;z0::AbstractVector = rand(length(y)), ϵ = 1e-2 ) -> z
Solves for a linear system of equations of the form `(S+D)z = y` iteratively. Here `S` is a sparse matrix, `D` is a dense matrix with rank 1 (PseudoDenseMatrix).
Not guaranteed to converge.
"""
function iterative_linsolve(S::SparseMatrixCSC,D::PseudoDenseMatrix,y::AbstractVector;z0::AbstractVector = rand(length(y)), ϵ = 1e-2,maxiter=1000)
	
	norm_diff = Inf
	z = z0
	i = 0
	
	while norm_diff > ϵ
		z = y + S*z - D*z
		norm_diff = norm(z-z0)
		
		norm_diff/length(z) > 1.0  && error("Norm of z-z0/length(z) is bigger than 1!")
		#@show norm_diff/length(z)
		
		i+=1
		i > maxiter && error("Maximum number of iterations exceeded! Your system converges slowly/might not converge at all. Try setting `maxiter` kwarg to bigger than 200, or `ϵ` to a higher value!")
		
		z0 = z
	end

	return z

end
