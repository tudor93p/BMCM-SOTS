module BBH 
#############################################################################

import LinearAlgebra 

import myLibs: Groups 
import myLibs.Parameters: UODict

import ..WLO#, ..MB  
import ..Helpers:Symmetries 
#import LinearAlgebra
#
#import ..WLO: kPOINT_START, NR_kPOINTS

const MODEL_MATRICES_PAULI::Vector{String} = ["30","21","22","23","10"]  

#const delta::Float64 = 1e-8

import ..MB: polariz_fromSubspaces1, add_nupm


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

usedkeys::Vector{Symbol} = [:BBHtheta]

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function H(k::AbstractVector{<:Real}, 
					 (lx, ly, gx, gy, delta)::AbstractVector{<:Real},
					 aux::Real=0
					 )::Matrix{ComplexF64}

return 

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


function BBHparams((th,)::AbstractVector{<:Real})::Vector{Float64} 

	BBHparams(th)

end  


function BBHparams(th::Real)::Vector{Float64} 

	BBHtheta = [-1, -0.04, 1.04, 2.04, 3, 3.96, 5, 6.04]/4

	@assert any(≈(th), pi*BBHtheta) th 

	[1, 1, 1 + 0.5*cos(th), 1+0.5*sin(th), 0]#1e-9]

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


function get_pertHdata(params::AbstractVector{<:Real}; kwargs...
							 )
	
	(BBHparams(params),)

end 

function get_pertHdata(params::AbstractVector{<:Real}, h::Function; kwargs...
							 )
	
	(h, get_pertHdata(params)...)

end 

function get_pertHdata(params::AbstractVector{<:Real}, p::AbstractArray; kwargs...
							 )
	
	(p, get_pertHdata(params)...)

end 

function get_pertHdata(params::AbstractVector{<:Real}, h::Function, p::AbstractArray; kwargs...
							 )
	
	(p, get_pertHdata(params, h)...)

end 

function get_pertHdata(params::AbstractVector{<:Real}, p::AbstractArray, s::Real;
							 atol::Float64=1e-12
							 )

	isapprox(s,0,atol=atol) && return get_pertHdata(params)  

	isapprox(s,1,atol=atol) && return get_pertHdata(params, p)

	return (get_pertHdata(params, p)..., s)

end 


function get_pertHdata(params::AbstractVector{<:Real}, h::Function, p::AbstractArray, s::Real;
							 atol::Float64=1e-12
							 )

	isapprox(s,0,atol=atol) && return get_pertHdata(params, h) 

	isapprox(s,1,atol=atol) && return get_pertHdata(params, h, p)
	
	return (get_pertHdata(params, h, p)..., s) 

end 


function symmetrize_on_mesh!(args...; kwargs...)
#NOT IMPLEMENTED 
end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_psiH(params::AbstractVector{<:Real}, n::Int, k0::Real, args...;
									atol::Float64=1e-12
								 )::Array{ComplexF64,4}

	WLO.psiH_on_mesh(n, k0, get_pertHdata(params, H, args...; atol=atol)...)

end 

function get_enpsiH(params::AbstractVector{<:Real}, n::Int, k0::Real, args...;
									atol::Float64=1e-12
									)::Vector{Array}

	WLO.enpsiH_on_mesh(n, k0, get_pertHdata(params, H, args...; atol=atol)...)

end 

function get_psiH(P::UODict, args...)::Array{ComplexF64,4}

	get_psiH(params_fromP(P), args...)

end  
function params_fromP(P::UODict)::Vector{Float64}

	[P[:BBHtheta]*pi, ]

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

R2repr = -im*Groups.GammaMatrix(0,2)

function operR2(h::AbstractMatrix)::Matrix 

	Symmetries.UdAV(R2repr, h)

end 

function operR2((V,U)::Tuple{AbstractMatrix,AbstractMatrix})::Matrix 

	Symmetries.UdAV(V, R2repr, U)

end 

function operR2(k::AbstractVector{<:Real})::Vector{Float64}

	Symmetries.inversion(k)

end 


function fij_operR2(n::Int, k0::Real)::Function # = Symmetries.inversion   

	ind_minusk = WLO.ind_minusk(n, k0)

	fij_(k::Vararg{Int,2})::NTuple{2,Int} = fij_(k) 

	fij_(k::NTuple{2,Int})::NTuple{2,Int} = Symmetries.inversion(k, ind_minusk)

end 



function uniqueInds_operR2(n::Int, k0::Real
																)::Tuple{<:AbstractVector{Int}, 
																				 <:AbstractVector{Int}
																				}

	(WLO.uniqueInds_kMirror(n,k0), Base.OneTo(n-1))

end  


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function get_perturb_on_mesh(
										 symms::AbstractString, 
										 n::Int, k0::Real,
										 )::Array{ComplexF64,4}

	@assert symms=="R2" 

	pert = rand(ComplexF64,4,4,n-1,n-1)
	
	WLO.store_on_mesh!!(Symmetries.symmetrize_HC!, pert)
		
	WLO.symmetrize_on_mesh!(pert,
													uniqueInds_operR2(n,k0),
													(operR2, fij_operR2(n,k0))
														)

#	MB.symmetrize_on_mesh!(pert, args...)

	WLO.store_on_mesh!!(Symmetries.normalize!, pert)
	
	return  pert 
end   




#############################################################################
end # module BBH
