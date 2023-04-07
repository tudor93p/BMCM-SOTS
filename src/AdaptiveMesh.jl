module AdaptiveMesh 
#############################################################################
#
# References
#	
#	"Numerical Recipes", Press et al.: 17.2 Adaptive Runge-Kutta p.934
#
# "Convergence of a step-doubling ...", B. Ayati, T.Dupont, 
# 																	Math. Comp. 74 (2005), 1053-1065#
#
#############################################################################




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function scale(y::Real, atol::Real, rtol::Real=atol)::Float64

	atol + y * rtol 

end 

function scale(y::Union{AbstractVector{<:Real},Tuple{Vararg{<:Real}}},
							 tols::Real...):Float64

	scale(maximum(abs,y), tols...)

end 




function timestep(timevalues::AbstractVector{Float64}, n::Int)::Float64
	
	timevalues[n]-timevalues[n-1] 

end 

function timestep(timevalues::AbstractVector{Float64}, n::Int,
								 factor::Real)::Float64

	timevalues(timevalues, n)*factor

end 



function new_timestep_failure(
					history::AbstractVector,
					 n::Int,
					 recommended_dt::Real,
					 desired_ts::AbstractVector{<:Real},
					 )

	foo(tjimestep(timevalues, n, 0.5))

end 


function bar!(
history::AbstractVector{Int},
timevalues::AbstractVector{Float64},
n::Int,
err::Real, tol::Real,
)

	if err>tol 

		history[n] -= 1 

		return (false,new_timestep_failure())

	else 

		history[n] 

		return (true, new_timestep_success())


	end 


end 







function step_selection( err::Real, tol::Real )

	err>tol && return (false,dt/2)

	err>3tol/4

end 










































































































































































































































































































#############################################################################
end # module AdaptiveMesh
