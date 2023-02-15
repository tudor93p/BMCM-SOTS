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

	ceil(max(0,(min(a,b)-j)/T))-ceil(max(0,(j-max(a,b))/T))

end 

function bring_periodic_to_interval(j::T1, a::Real, b::Real, T::T2
																		)::promote_type(T1,T2) where {
																												T1<:Real,T2<:Real}

	j + T*periods_outside_interval(j, a, b, T)

end  




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






















#############################################################################
end # module PeriodicFuns

