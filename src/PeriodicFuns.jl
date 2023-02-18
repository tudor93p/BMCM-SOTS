module PeriodicFuns 
#############################################################################

import myLibs:Utils 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function closest_periodic_a(a::Real,
													b::Real,
													T::Real,
													nmax::Int=1,
													)::Float64 

#	a + T*argmin((n::Int)::Float64->abs(a + n*T - b), -nmax:nmax)
	a_::Float64 = copy(a)

	for n in -nmax:nmax

		if abs(a + n*T - b) < abs(a_ - b)

			a_ = a + n*T 

		end 

	end 

	return a_

end  

function closest_periodic_a(a::Real,
													B::AbstractArray{<:Real},
													T::Real,
													nmax::Int=1
													)::Float64 

	out = Float64[a, abs(a-B[1])]
	
	current = Float64[a+nmax*T, 0]

	for n in -nmax:nmax

		for b in B 
	
			current[2] = abs(current[1]-b)
			
			if current[2] < out[2]
				
				out .= current
	
			end 
	
		end 

		current[1] -= T  

	end 

	return out[1] 

end 


function closest_periodic_b(a::Real,
														B::AbstractArray{<:Real},
														args...
													)::Float64

	B[argmin(Utils.dist_periodic(a, B, args...))]

end 
#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function periods_outside_interval(j::Real, a::Real, b::Real, T::Real)::Int 

#	ceil(max(0,(a-j)/T))-ceil(max(0,(j-b)/T))

	j<=a && return  ceil((a-j)/T) 
	j>=b && return -ceil((j-b)/T)

end 

function bring_periodic_to_interval(j::T1, a::Real, b::Real, T::T2
																		)::promote_type(T1,T2) where {
																												T1<:Real,T2<:Real}
	j + T*periods_outside_interval(j, a, b, T)

end  

#
#function bring_periodic_to_interval!(container::AbstractArray{<:Real}, 
#																		 j::Int, 
#																		 a::Real, b::Real, T::Real
#																		)::Nothing 
#
#	container[j] += T*periods_outside_interval(container[j], a, b, T)
#
#	return 
#
#end  



#function dist_modulo(x::Real,m::Real=1)::Real 
#
#	min(abs(x),abs(m-abs(x)))
#
#end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function reduce_index(i::Int, period::Int)::Int 

	i>period && return reduce_index(i-period, period)

	i<1 && return reduce_index(i+period, period)

	return i

end 


function cycseldim(A::AbstractArray{T,N}, d::Int, i::Int
									 )::AbstractArray{T,N-1} where {N,T<:Number}

	selectdim(A, d, reduce_index(i, size(A,d)))

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function dist_periodic_!(D::AbstractArray{Float64},
												 i::Int,
												 a::Real,
												 b::Real,
												 T::Real,
												 )::Nothing 

	setindex!(D, bring_periodic_to_interval(a-b,0,T,T), i)

	setindex!(D, min(D[i], T-D[i]), i)

	return 

end  


function reduce_dist_periodic(op::Function,
															A::Union{<:Real, AbstractArray{<:Real}},
															B::Union{<:Real, AbstractArray{<:Real}},
															args...; init=nothing)

	isempty(A) && return init 
	isempty(B) && return init 

	hasmethod(op, (Float64,)) || return init  

	(isa(A,AbstractArray) & isa(B,AbstractArray)) && @assert size(A)==size(B)
	
	D = Vector{Float64}(undef, 1)

	dist_periodic_!(D, 1, get(A,1,A), get(B,1,B), args...)

	out = [op(only(D))]

	length(A)==1 && length(B)==1 && return only(out)

	for i=2:max(length(A),length(B))

		dist_periodic_!(D, 1, get(A,i,A), get(B,i,B), args...)

		setindex!(out, op(only(out), only(D)), 1) 

	end 

	return only(out)

end 









#############################################################################
end # module PeriodicFuns

