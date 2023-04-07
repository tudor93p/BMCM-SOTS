module CalcWLOadapt
#############################################################################

import Random, LinearAlgebra


using OrderedCollections: OrderedDict 

import myLibs: ReadWrite, Utils

import myLibs.Parameters: UODict 

import ..FILE_STORE_METHOD


import ..MB; MODEL=MB
##import ..BBH; MODEL=BBH 


import ..WLO
import ..WLO: nr_kPoints, kPoint_start  

import ..Helpers: Symmetries 


import ..CalcWLO
import ..CalcWLO: all_symms_preserved, preserved_symmetries,
									 perturb_strength, 
									 islegend, fn_legend, xxlabel,
									 get_target, FoundFiles, Read,
									 cp_wcc_to_results,
									 get_data_args,
									 set_results!


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

Dependencies = [MODEL,WLO,CalcWLO]

usedkeys::Vector{Symbol} = [:kMesh_type]


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function init_results(n::Int, k0::Real,
											obs::AbstractVector{<:AbstractString}
										 )::Dict{String,Any}

	results = CalcWLO.init_results(n, obs) 

# "ks" => WLO.get_kij(n, k0; restricted=false)(1:n)

	@warn "'ks' not yet in 'results'"

		#@assert Set(keys(results))==Set(xxlabel())

	return results

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


function symmetrize_and_normalize!(pert::AbstractArray{ComplexF64,4}, 
																	 symms::AbstractString,
																	 args...)::AbstractArray{ComplexF64,4}

	WLO.store_on_mesh!!(Symmetries.symmetrize_HC!, pert)

	if !isempty(MODEL.sepSymmString(symms))

		error("Not ready to symmetrize on non-uniform mesh")

	end 
#	MODEL.symmetrize_on_mesh!(pert, symms, args...)

	WLO.store_on_mesh!!(Symmetries.normalize!, pert)

	return pert 

end 



function get_perturb_on_mesh_(
										 symms::AbstractString, 
										 )::Array{ComplexF64,4}

	symmetrize_and_normalize!(rand(ComplexF64,4,4,n-1,n-1), symms, n, k0) 

end   

function get_perturb_on_mesh_(
										 symms::AbstractString, 
										 n::Int, k0::Real,
										 seed::Int 
										 )::Array{ComplexF64,4}


	Random.seed!(seed) 

	return get_perturb_on_mesh_(symms, n, k0)
 
end 

function get_perturb_on_mesh(symms::AbstractString, n::Int,
														 args...)::Array{ComplexF64,4}
	
	all_symms_preserved(symms) && return zeros(ComplexF64,4,4,n-1,n-1)

	return get_perturb_on_mesh_(symms, n, args...)

end 

function get_perturb_on_mesh(
											P::UODict, args...
										 )::Array{ComplexF64,4}

	get_perturb_on_mesh(preserved_symmetries(P), 
											nr_kPoints(P),
											kPoint_start(P),
											args... 
											) 

end 













#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function get_data(psiH::AbstractArray{ComplexF64,4}, 
									results::AbstractDict;
									parallel::Bool=false
									)::Vector

	error("'get_data' not implemented")

	#(parallel ? pmap : map)(Base.splat(WLO.get_wlo_data_mesh), 
	#												get_data_args(psiH, results))

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


#===========================================================================#
#
# Compute 
#
#---------------------------------------------------------------------------#




function Compute(P::UODict; get_fname=nothing, target=nothing,
								 kwargs...)::Dict  

	Compute_(P, target, get_fname; kwargs...)

end   




function Compute_(P::UODict, target, get_fname::Function;
									kwargs...)::Dict{String,Any}

	obs = get_target(target; kwargs...) 
	
	results = Compute_(P, obs; kwargs...)

	ReadWrite.Write_PhysObs(get_fname(P), FILE_STORE_METHOD, results)

	return isnothing(target) ? results : Utils.dict_keepkeys(results, obs) 

end 





function Compute_(P::UODict, target, get_fname::Nothing=nothing; 
										kwargs...)::Dict{String,Any}

	nk = nr_kPoints(P)
	k0 = kPoint_start(P)
	symms = preserved_symmetries(P)
	strength = perturb_strength(P)

	results = init_results(nk, k0, get_target(target; kwargs...))

	@show results


	perturb1 = get_perturb_on_mesh(P)

	@show LinearAlgebra.norm(perturb1)

	nk_unif = max(5,div(nk,3))

	psi = MODEL.get_psiH(P, nk, k0, perturb1, strength)



#	set_results!(results, nk, k0, get_data(psi, results))

	return results

end 







































































































#############################################################################
end # module CalcWLO

