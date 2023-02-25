module BBH 
#############################################################################

import LinearAlgebra 

import myLibs: Groups
import ..WLO#, ..MB 
#import LinearAlgebra
#
#import ..WLO: kPOINT_START, NR_kPOINTS

const MODEL_MATRICES_PAULI::Vector{String} = ["30","21","22","23","10"]  

const delta::Float64 = 1e-8

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function H(k::AbstractVector{<:Real}, 
					 (lx, ly, gx, gy, delta)::AbstractVector{<:Real},
					 aux::Real=0
					 )::Matrix{ComplexF64}

	sx,cx = sincos(k[1])
	sy,cy = sincos(k[2])


	P = [delta, ly*sy, gy+ly*cy, lx*sx, gx+lx*cx] 

	P[2:4] *= -1  # Kron def without minus 

	out = zeros(ComplexF64, 4, 4) 

	G = zeros(ComplexF64, 4, 4) 

	for (p,ij) in zip(P,MODEL_MATRICES_PAULI)

		LinearAlgebra.axpy!(p, Groups.GammaMatrix!(G, ij), out)

	end 

	return out

end 

function H(k::AbstractVector{<:Real},
					 bs::AbstractVector{<:Real},
					 perturb::AbstractMatrix{<:Number},
					 )::Matrix{ComplexF64}
	
	H(k, bs) .+= perturb 

end  


function H(k::AbstractVector{<:Real},
					 bs::AbstractVector{<:Real},
					 pstrength::Real, 
					 perturb0::AbstractMatrix{<:Number}, 
					 )::Matrix{ComplexF64}

	LinearAlgebra.axpy!(pstrength, perturb0, H(k,bs))

end  


function H(ij::NTuple{2,Int},
						kpoint::Union{Function,<:AbstractArray{<:Real,4}},
						perturb::Union{Function,<:AbstractArray{<:Number,4}},
						args...
					 )::Matrix{ComplexF64}

	H(WLO.get_item_ij(ij, kpoint), args..., WLO.get_item_ij(ij, perturb))

end 

function H(ij::NTuple{2,Int},
						kpoint::Union{Function,<:AbstractArray{<:Real,4}},
						bs::AbstractVector{<:Real},
					 )::Matrix{ComplexF64}

	H(WLO.get_item_ij(ij, kpoint), bs)

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function BBHparams(th::Real)::Vector{Float64} 

	[1, 1, 1 + 0.5*cos(th), 1+0.5*sin(th), 1e-8]

end  



function BBHoutcomes(th::Real)::Vector{Float64}#String

	lx, ly, gx, gy, = BBHparams(th)

	return [abs(gx/lx)<1, abs(gy/ly)<1] * 0.5 

	#	Q = [(1,1),(1,0),(0,1),(0,0),(0,1), (0,1)][i]
end 

function BBHoutcomes_str(th::Real)::String

	string(
				rpad(string(round(th/pi,digits=3)),5,'0'),
				"\$\\pi\$: ",
				"(",
				join([q≈0.5 ? "\$\\pm 0.5\$" : 0 for q=BBHoutcomes(th) 
							],", "),
				")"
				)
end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_pertHdata(BBHtheta::Real; kwargs...
							 )
	
	(BBHparams(BBHtheta),)

end 

function get_pertHdata(BBHtheta::Real, h::Function; kwargs...
							 )
	
	(h, get_pertHdata(BBHtheta)...)

end 

function get_pertHdata(BBHtheta::Real, p::AbstractArray; kwargs...
							 )
	
	(p, get_pertHdata(BBHtheta)...)

end 

function get_pertHdata(BBHtheta::Real, h::Function, p::AbstractArray; kwargs...
							 )
	
	(p, get_pertHdata(BBHtheta, h)...)

end 

function get_pertHdata(BBHtheta::Real, p::AbstractArray, s::Real;
							 atol::Float64=1e-12
							 )

	isapprox(s,0,atol=atol) && return get_pertHdata(BBHtheta)  

	isapprox(s,1,atol=atol) && return get_pertHdata(BBHtheta, p)

	return (get_pertHdata(BBHtheta, p)..., s)

end 


function get_pertHdata(BBHtheta::Real, h::Function, p::AbstractArray, s::Real;
							 atol::Float64=1e-12
							 )

	isapprox(s,0,atol=atol) && return get_pertHdata(BBHtheta, h) 

	isapprox(s,1,atol=atol) && return get_pertHdata(BBHtheta, h, p)
	
	return (get_pertHdata(BBHtheta, h, p)..., s) 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_psiH(BBHtheta::Real, n::Int, k0::Real, args...;
									atol::Float64=1e-12
								 )::Array{ComplexF64,4}

	WLO.psiH_on_mesh(n, k0, get_pertHdata(BBHtheta, H, args...; atol=atol)...)

end 






#############################################################################
end # module BBH
