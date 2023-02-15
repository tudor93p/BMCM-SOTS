module Algebra 
#############################################################################



#===========================================================================#
#
# more generically implemented Modified Gram-Schmidt. Update myLibs.Algebra
#
#---------------------------------------------------------------------------#
#
#function GS_get_ind_matr(X::AbstractMatrix, i::Int; dim::Int, kw...
#								 )::AbstractVector
#
#		selectdim(X, dim, i)
#
#end 
#
#
#function GramSchmidt!(V::AbstractMatrix, iter_axis::Int=2;
#										 kwargs...)::Nothing
#
#	GramSchmidt_!(V, V, size(V,iter_axis), GS_get_ind_matr; 
#								dim=iter_axis, kwargs...)
#
#end 
#
#
#function GramSchmidt!(V::AbstractVector{<:AbstractVector};
#										 kwargs...
#										 )::Nothing 
#
#	GramSchmidt_!(V, V, length(V), getindex; kwargs...)
#
#end 
#
#
#
#
#function GramSchmidt(V::AbstractVector{<:AbstractVector}; kwargs...
#										 )::Vector{Vector{<:Number}}
#
#	isempty(V) && return Vector{Vector{Float64}}(undef, 0)
#
#	isempty(V[1]) && return [Vector{Float64}(undef, 0) for v in V]
#
#
#	U = [zeros(promote_type(Float64, typeof(V[1][1])), length(V[1])) for v in V]
#
#	GramSchmidt_!(U, V, length(V), getindex; kwargs...)
#
#	return U
#
#end 
#
#
#function GramSchmidt(V::AbstractMatrix, iter_axis::Int=2;
#										 kwargs...)::Matrix{<:Number} 
#
#	isempty(V) && return Matrix{Float64}(undef, 0, 0)
#
#
#	U = zeros(promote_type(Float64, typeof(V[1])), size(V))
#
#	GramSchmidt_!(U, V, size(V,iter_axis), GS_get_ind_matr; 
#								dim=iter_axis, kwargs...)
#
#	return U 
#
#end 
#
#
#
#function GramSchmidt_!(out_U::Tu, inp_V::Tv,
#											 jmax::Int,
#											 get_ind::Function;
#											 normalize::Bool=true, tol::Float64=1e-10,
#											 kwargs...
#											 )::Nothing where {
#																				 T<:Union{AbstractMatrix, 
#																									AbstractVector{<:AbstractVector}},
#																				 Tu<:T, Tv<:T} 
#
#
#	jmax < 1 && return 
#
#
#	u = zero(get_ind(out_U,1; kwargs...))
#
#	@assert jmax<=length(u) "Too many vectors"
#
#	norm = zeros(1)
#
#	for j=1:jmax 
#
#		setindex!(u, get_ind(inp_V, j; kwargs...), :)
#
#		for i in 1:j-1
#
#			u -= LinearAlgebra.dot(get_ind(out_U, i; kwargs...), u) * get_ind(out_U, i; kwargs...)
#
#		end 
#
#
#		if normalize 
#		
#			norm[1] = LinearAlgebra.norm(u)
#
#			@assert norm[1]>tol "Invalid input matrix for Gram-Schmidt!"
#			
#			u ./= norm[1]
#
#		end 
#
#		setindex!(get_ind(out_U, j; kwargs...), u, :)
#
#	end
#
#  return
#
#end
#
#

#############################################################################


end 
