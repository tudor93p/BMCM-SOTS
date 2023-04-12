module AdaptiveMesh 
#############################################################################

import myLibs: Utils 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function f_line(p::AbstractVector{<:Real}
								)::Function 

	f(x::Real)::Float64 = p[1]*x+p[2] 

end 

function line_par_from_vals(X::AbstractVector{<:Real},
														Y::AbstractVector{<:Real},
														)::Vector{Float64}

#	@assert all(>(0),X)&&all(>(0),Y) "Non-negative values expected"

	@assert !(X[1]≈X[2]) "X values to close from eachother"

	Y[1]≈Y[2] && return zeros(2) 

	@assert !xor(X[1]<X[2],Y[1]<Y[2]) "The function must increase"

	return inv(hcat(X,ones(2)))*Y

end 

function f_square(p::AbstractVector{<:Real}
								)::Function 

	abs2∘f_line(p)

end 

function f_sin(p::AbstractVector{<:Real})::Function 

	l = f_line(view(p,1:2))

	f(x::Real)::Float64 = p[3]*sinpi(l(x))

	return f 

end 

#function f_cos(p::AbstractVector{<:Real})::Function 
#
#	l = f_line(view(p,1:2))
#
#	f(x::Real)::Float64 = p[3]*(1-cospi(l(x)))
#
#	return f 
#
#end 


function f_exp(p::AbstractVector{<:Real})::Function 

	exp∘f_line(p)

end 

function exp_par_from_vals(X::AbstractVector{<:Real},
														Y::AbstractVector{<:Real},
														)::Vector{Float64}

	line_par_from_vals(X,log.(Y))

end 
function sin_par_from_vals(X::AbstractVector{<:Real},
														Y::AbstractVector{<:Real},
														)::Vector{Float64}

	ym,yM = sort(Y)

	return vcat(line_par_from_vals(X, [asin(ym/yM)/pi,0.5]), yM)

end  
#function cos_par_from_vals(X::AbstractVector{<:Real},
#														Y::AbstractVector{<:Real},
#														)::Vector{Float64}
#
#	ym,yM = sort(Y)
#
#	return vcat(line_par_from_vals(X, [1/pi-acos(ym/yM)/pi,0.5]), yM)
#
#end  


function expminv_par_from_vals(X::AbstractVector{<:Real},
														Y::AbstractVector{<:Real},
														)::Vector{Float64}

	line_par_from_vals(X,-inv.(log.(Y)))

end 

function square_par_from_vals(X::AbstractVector{<:Real},
														Y::AbstractVector{<:Real},
														)::Vector{Float64}

	line_par_from_vals(X,sqrt.(Y))


end 

function f_expminv(p::AbstractVector{<:Real})::Function 

	∘(exp, -, inv, f_line(p))

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





bound_rescale_kStep(dk::Real)::Float64 = dk 



function bound_rescale_kStep(dk::Real,
											(dk_min,dk_max)::AbstractVector{<:Real},
											alpha::Real=1
											)::Float64

	max(min(alpha*dk,dk_max),dk_min)

end 


#function smoothmax(a::Real,b::Real)::Float64 
#
#	@assert !xor(a<0, b<0)
#	
#	
#	v = [a,b] 
#
##	sign(argmax(abs, v))
#
#	return sign(a)*LinearAlgebra.norm(v,6)
#
#end 
#
#function smoothmin(a::Real,b::Real)::Float64 
#
#	a+b-smoothmax(a,b) 
#
#end 
#
#softbound(dk::Real, dk_min::Real, dk_max::Real)::Float64 = smoothmin(smoothmax(dk,dk_min),dk_max)
#
#hardbound(dk::Real, dk_min::Real, dk_max::Real)::Float64 = max(min(dk,dk_max),dk_min)
#
#
#
#
#function softBound_rescale_kStep(dk::Real,
#											bounds::AbstractVector{<:Real},
#											alpha::Real=1
#											)::Float64
#
#	softbound(dk*alpha, bounds...)
#
#end 
#
#softBound_rescale_kStep(dk::Real,)::Float64 = dk 
#
#
#function softBound_rescale_kStep(get_dk::Function,
#														 br_args...
#											)::Function 
#
#	function get_dk_(args...)::Float64
#	
#		softBound_rescale_kStep(get_dk(args...), br_args...)
#
#	end 
#
#end 
#softBound_rescale_kStep = bound_rescale_kStep


function bound_rescale_kStep(get_dk::Function,
														 br_args...
											)::Function 

	function get_dk_(args...)::Float64
		
		bound_rescale_kStep(get_dk(args...), br_args...)

	end 

end 
	

function verify_dk_bounds(
											(dk_min,dk_max)::AbstractVector{<:Real},
											N::Int
											)::Vector{Float64}

	K = [dk_min, dk_max]

	if dk_min*(N-1)>2pi 
		@warn "dk_min too large. Using 2pi/(N-1)"

		K[1] = 2pi/(N-1)

	end 

	if dk_max*(N-1)<2pi 

		@warn "dk_max too small. Using 2pi/(N-1)"

		K[2] = 2pi/(N-1)
	end 

	return K

end  

#
#function find_largest_step(ks::AbstractVector{<:Real},
#													 )::Tuple{Float64,Int}
#
#	i0::Int = 0 
#	dk0::Float64 = 0.0
#	dk::Float64 = 0.0
#
#	for i = 2:length(ks)
#
#		dk = Utils.dist_periodic(ks[i-1],ks[i],2pi)
#
#		if dk > dk0
#
#			i0 = i
#			dk0 = dk 
#
#		end
#	end 
#
#	return dk0,i0 
#
#end 
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#






































#############################################################################
#
# References for solving diff. eq.
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

	foo(timestep(timevalues, n, 0.5))

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

function baz(f::Function, 
						 data,
						 x0::Float64
						 )

	f(data, x)

end 


function step_selection( err::Real, tol::Real )

	err>tol && return (false,dt/2)

	err>3tol/4

end 










































































































































































































































































































#############################################################################
end # module AdaptiveMesh
