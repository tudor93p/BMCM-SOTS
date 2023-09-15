module ChiralPwave 

import LinearAlgebra 
#, Statistics 
import SharedArrays: SharedArray
#
import ..WLO 
#import ..Helpers:Symmetries 
##import PyPlot 
import myLibs: Groups, Utils

import myLibs.Parameters: UODict

#import DelimitedFiles 

const quantized_wcc2_values = [0,0.5]

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

needs_config(P::UODict)::Bool = P[:pairSymm]=="chiral"


usedkeys()::Vector{Symbol} = [:Delta0,
														:pairSymm,
														:pairSymmConfig,
														]

function usedkeys(P)::Vector{Symbol}
	
	needs_config(P) ? usedkeys() : setdiff!(usedkeys(),[:pairSymmConfig])

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function dispersion((kx,ky)::AbstractVector{<:Real})::Float64

	2 - 2cos(kx) - 2cos(ky)

end 


function dvector_chiral((kx,ky)::AbstractVector{<:Real},
												(D0, CN,)::AbstractVector{<:Real},
								 )::Vector{ComplexF64}

	ComplexF64[0, 0, D0*(sin(kx) + im*CN*sin(ky))] 

end 


function psi_extswave((kx,ky)::AbstractVector{<:Real},
												(D0, )::AbstractVector{<:Real},
								 )::ComplexF64

	D0*cos(kx)+D0*cos(ky)

end 



function dvector_helical((kx,ky)::AbstractVector{<:Real},
												 (D0,)::AbstractVector{<:Real},
								 )::Vector{ComplexF64}

	ComplexF64[D0*sin(ky), -D0*sin(kx), 0] 

end 


function gap_function!(Delta::T,
											 d::AbstractVector{<:Number}
											 )::T where T<:AbstractMatrix{ComplexF64}

	@assert size(Delta)==(2,2)

	Delta .= 0

	for (i,di) in enumerate(d) 

		LinearAlgebra.axpy!(im*di, Groups.PMprod(i,2), Delta)

	end 

	return Delta 

end 

function gap_function!(Delta::T,
											 psi::Number
											 )::T where T<:AbstractMatrix{ComplexF64}

	@assert size(Delta)==(2,2)

	Delta .= Groups.PauliMatrix(2)

	Delta .*= im*psi 

	return Delta 

end 





function separate_H_diagblocks(h::AbstractMatrix{T},
															 I1::AbstractVector{Int};
															 atol::Float64=1e-6,
															 )::NTuple{2,AbstractMatrix{T}
																				 } where T<:Number 

	I2 = setdiff(axes(h,1), I1)

	for i1=I1, i2=I2 

		isapprox(h[i1,i2], 0, atol=atol) && continue 
		isapprox(h[i2,i1], 0, atol=atol) && continue 

		@warn "Blocks $I1 $I2 coupled!"

		break 

	end 

	return view(h, I1, I1), view(h, I2, I2) 

end 


function H(k::AbstractVector{<:Real}, 
					 gap_args::AbstractVector{<:Real},
					 psi_or_d::Function,
					)::Matrix{ComplexF64}

	delta = psi_or_d(k, gap_args)


	h = zeros(ComplexF64, 4, 4)

	h[1,1] = h[2,2] = dispersion(k) 
	h[3,3] = h[4,4] = -h[1,1]

	view(h, 3:4, 1:2) .= gap_function!(view(h,1:2,3:4), delta)'

	return h 

#	h1, h2 = separate_H_diagblocks(h, [1,4])
#
#	@assert isapprox(h1,h2,atol=1e-6) 
#
#	return h1

end 



#function H(ij::NTuple{2,Int},
#						kpoint::Union{Function,<:AbstractArray{<:Real,4}},
#						params::AbstractVector{<:Real},
#					 )::Matrix{ComplexF64}
#
#	H(WLO.get_item_ij(ij, kpoint), params)
#
#end 
#


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function params_fromP(P::UODict)::Tuple{Vector{Float64},Function}

	@show P 

	gap_magnitude = P[:Delta0]

	if P[:pairSymm]=="chiral"

		chern_number = [-1,1][P[:pairSymmConfig]]

		return ([gap_magnitude, chern_number], dvector_chiral)

	elseif P[:pairSymm]=="helical"

		return ([gap_magnitude], dvector_helical) 

	elseif P[:pairSymm]=="ext.swave"

		return ([gap_magnitude], psi_extswave)

	else 

		error("Not implemented")

	end 

end 






function get_psiH(P::UODict,
									n::Int, 
									k::Union{<:Real,<:AbstractVector{<:AbstractVector{<:Real}}
														},
									args...; # perturbation ignored 
									kwargs...
								 )::Union{Array{ComplexF64,4},
													SharedArray{ComplexF64,4}
													}


	WLO.psiH_on_mesh(n, k, H, params_fromP(P)...; kwargs...)

end  







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function add_nupm(storage::AbstractVector{<:AbstractArray},
#									dir::Int; kwargs...
#								 )::Vector{Float64}
#
#	out = zeros(size(storage[1],2))
#
#	for s in storage 
#
#		@assert ndims(s)==3
#
#		out .+= WLO.check_nu_k_dep(selectdim(s,1,1), dir; kwargs...)[2]
#
#	end 
#
#	return out 
#		
#end    
#
#function polariz_fromSubspaces1(storage::AbstractVector{<:AbstractArray},
#																dir::Int; kwargs...)::Vector{Float64}
#
#	add_nupm(view(storage,[2,4]), dir; kwargs...)
#
#end  
#





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#






#############################################################################
end 
