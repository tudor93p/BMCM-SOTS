module Symmetries 
#############################################################################

import LinearAlgebra

import myLibs: Operators

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function UdAV(
							U::AbstractMatrix{Tu},
							A::AbstractMatrix{Ta},
							V::AbstractMatrix{Tv}=U,
							)::Matrix{promote_type(Tu,Ta,Tv)
												} where {T<:Number,Tu<:T,Ta<:T,Tv<:T}

	Operators.MatrixElem(U, A, V, 1, :) # U'*A*U'

end  


#function apply_unit_symm(repr::AbstractMatrix{Tu},
#							A::AbstractMatrix{Ta}
#							)::Matrix{promote_type(Tu,Ta)} where {Tu<:Number,Ta<:Number}
#
#
#	UdAU(repr,A) 
#
#end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function mirror(k::AbstractVector{T}, d::Int)::Vector{T} where T<:Real 

	setindex!(copy(k), -k[d], d)

end 

function mirror(k::NTuple{N,Int}, d::Int, rev::Function)::NTuple{N,Int} where N

	Tuple(i==d ? rev(k[i]) : k[i] for i=1:N)

end 

function inversion(k::AbstractVector{T})::Vector{T} where T<:Real 

	inversion!(copy(k))

end  

function inversion(k::NTuple{N,Int}, rev::Function)::NTuple{N,Int} where N 

	map(rev, k)

end   

function inversion!(k::AbstractVector{T})::Vector{T} where T<:Real 

	k .*= -1 

end  



function rotationC2(k::AbstractVector{T}, d::Int)::Vector{T} where T<:Real 

	[d==i ? -ki : ki for (i,ki)=enumerate(k)]

end 

function rotationC2(k::NTuple{N,Int}, d::Int, rev::Function
										)::NTuple{N,Int} where N 

	Tuple(d==i ? rev(ki) : ki for (i,ki)=enumerate(k))

end 









#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function symmetrize(A::AbstractMatrix, args...)::Matrix

	B = copy!(zeros(ComplexF64,size(A)),A)

	symmetrize!(B, args...)

	return B 

end 



########

#function symmetrize!(aij::AbstractMatrix,
#										 op::Function,
#										 Amesh::AbstractArray{ComplexF64,4},
#										 ij::NTuple{2,Int},
#										 ind_minusk::Function,
#							 )::Nothing 
#
#	A .+= op(select_mesh_point(Amesh, op(i,j,ind_minusk)))
#
#end 

#function symmetrize!(A::AbstractMatrix, op::Function, 
#										 )::Nothing
#
#	A .+= op(A)
#
#	return
#
#end   
##########



function symmetrize!(A::AbstractMatrix, op::Function, 
										 B::AbstractMatrix=A
										 )::Nothing

	A .+= op(B)

	return

end   


#mutual symmetrization, first two entries modified in-place   

function symmetrize!!(A::AbstractMatrix, B::AbstractMatrix, op::Function
										 )::Nothing

	A .+= op(B) 

	B .= op(A) 

	return 

end    



function symmetrize_HC!(A::AbstractMatrix, args...)::Nothing

	symmetrize!(A, adjoint, args...)

end 

function normalize!(A::AbstractMatrix, args...)::Nothing

	LinearAlgebra.normalize!(A)

	return 

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function has_symm(op::Function, A::AbstractMatrix, B::AbstractMatrix=A;
									atol::Float64=1e-8)::Bool

	isapprox(A,op(B),atol=atol)

end   



























































#############################################################################
end # module Symmetries


