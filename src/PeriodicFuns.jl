module PeriodicFuns 
#############################################################################

import myLibs:Utils 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#
function closest_periodic_a2(a::Real,
													b::Real,
													T::Real,
													nmax::Int=1,
													)::Float64 

	a_::Float64 = copy(a)

	for n in -nmax:nmax

		if abs(a + n*T - b) < abs(a_ - b)

			a_ = a + n*T 

		end 

	end 

	return a_

end   

closest_periodic_b2(a,B,T) =B[argmin(Utils.dist_periodic(a,B,T))]

function closest_periodic_a2(a::Real,
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






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



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









#############################################################################
end # module PeriodicFuns

