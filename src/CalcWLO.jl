module CalcWLO 
#############################################################################

using Distributed 
import Random 

import myLibs: ReadWrite, Utils

import myLibs.Parameters: UODict 

import ..FILE_STORE_METHOD

import ..MB, ..WLO, ..Helpers

import ..MB: braiding_time
import ..Helpers: Symmetries 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


usedkeys = [
						:braiding_time,
						:nr_kPoints,
						:kPoint_start,
						:preserved_symmetries,
						:perturb_strength,
						]
									 

calc_observables = ["WannierGap","WannierBands1","WannierBands2"]


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function nr_kPoints(P::UODict)::Int 

	P[:nr_kPoints]

end  

function kPoint_start(P::UODict)::Float64

	pi*P[:kPoint_start]

end  
function preserved_symmetries(P::UODict)::String 

	P[:preserved_symmetries]

end 


function all_symms_preserved(P::UODict)::Bool

	all_symms_preserved(preserved_symmetries(P))

end 
function all_symms_preserved(s::AbstractString)::Bool

	s=="All"

end 

function perturb_strength(P::UODict)::Float64

	P[:perturb_strength]

end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#function xxlabel()::Vector{String}
#
#	["xs1","xlabels"]
#
#end  
#
#isauxfile = isxxlabel = in(xxlabel())
#
#function xxlabel(data::AbstractDict)::Tuple{<:AbstractMatrix,
#																						<:AbstractVector{<:AbstractString}
#																						}
#
#	Tuple(getindex(data,k) for k=xxlabel())
#
#end 


prep_obs(::Nothing=nothing)::Vector{String} = String[]

function prep_obs(obs::AbstractString)::Vector{String} 

	#obs in calc_observables ? vcat(obs,xxlabel()) : prep_obs() + legend? 

	obs in calc_observables ? vcat(obs) : prep_obs()

end 

function prep_obs(obs::AbstractVector{<:AbstractString})::Vector{String} 

	mapreduce(prep_obs, vcat, obs; init=prep_obs())

end 

function prep_obs(args...)::Vector{String} 
	
	mapreduce(prep_obs, vcat, args; init=prep_obs())

end 



get_target = prep_obs ∘ Helpers.f_get_target(:observables)	   




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function init_results(L::Int,
											obs::AbstractVector{<:AbstractString}
										 )::Dict{String,Any}


	if "WannierBands2" in obs 

		union!(obs, calc_observables)

	elseif "WannierBands1" in obs 

		union!(obs, ["WannierGap"])

	end 

	results = Dict{String,Any}()
	#"xs" => strengths, "xlabels"=>"kx ky etc."



	for w in ("WannierBands1","WannierBands2")

		in(w,obs) || continue 

		results[w] = zeros(L,2) 

	end 

	for w in ("WannierGap",)
		
		in(w,obs) || continue 

		results[w] = zeros(2) 

	end 

	return results

end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#
function symmetrize_on_mesh_one!(Aij::AbstractMatrix{ComplexF64},
									opij::NTuple{2,Int}, 
									seed::AbstractArray{<:Number,4},
									op::Function; kwargs...)::Nothing

	Symmetries.symmetrize!!(Aij, WLO.select_mesh_point(seed, opij), op)

end 


function symmetrize_on_mesh(seed::AbstractArray{ComplexF64,4},
														args...; kwargs...
														 )::Array{ComplexF64,4}

	A = copy(seed)

	symmetrize_on_mesh!(A, args...; kwargs...)

	return A

end  


function symmetrize_on_mesh!(
														 A::AbstractArray{ComplexF64,4},
														 opers::AbstractString,
														 n::Int, k0::Real,
														 )::Nothing 

	for op in MB.sepSymmString(opers)

		WLO.store_on_mesh!!(symmetrize_on_mesh_one!,
									MB.getOp_uniqueInds(op,n,k0),
									MB.getOp_fij(op,n,k0),
									A, A, 
									MB.getOpFun(op),
									) 
	end 

end  


function symmetrize_and_normalize!(pert::AbstractArray{ComplexF64,4}, 
																	 args...)::AbstractArray{ComplexF64,4}

	WLO.store_on_mesh!!(Symmetries.symmetrize_HC!, pert)

	symmetrize_on_mesh!(pert, args...)

	WLO.store_on_mesh!!(Symmetries.normalize!, pert)

	return pert 

end 



function get_perturb_on_mesh_(
										 symms::AbstractString, 
										 n::Int, k0::Real,
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
											) #.*= perturb_strength(P)

end 













#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_data_args(psiH::AbstractArray{<:Number,4},
											 results::AbstractDict
											 )::Vector{Tuple{Array,Bool,Int,Bool}}

	r = "WannierBands2" in keys(results)

	# only occupied subspace 
	return [(psiH, true, 1, r), #(psiH, false, 1, r), 
					(psiH, true, 2, r), #(psiH, false, 2, r), 
					]

end  

function get_data(psiH::AbstractArray{ComplexF64,4}, 
									results::AbstractDict;
									parallel::Bool=false
									)::Vector

	(parallel ? pmap : map)(Base.splat(WLO.get_wlo_data_mesh), 
													get_data_args(psiH, results))

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function check_nu_k_dep(nu::AbstractMatrix{<:Real},
												d::Int;
												atol::Float64=1e-10
												)

	for i=2:size(nu,d) 
	
		isapprox(selectdim(nu,d,i), selectdim(nu,d,1); atol=atol) && continue 
	
		@warn "nux depends on kx | nuy depends on ky?" 
	
		@show LinearAlgebra.norm(selectdim(nu,d,i)-selectdim(nu,d,1))

		break 
	
	end  

end 

function cp_wcc_to_results(results::AbstractDict{<:AbstractString,<:Any},
													K::AbstractString,
													wccs::Union{Tuple{<:T,<:T}, AbstractVector{<:T}
																		 } where T<:AbstractArray{<:Real,3},
													d::Int,
													n::Int
													)

	for (sector,nu_sector) in enumerate(wccs)

		check_nu_k_dep(selectdim(nu_sector,1,1), d) # one band 

		for k=1:n 

			setindex!(results[K],
								only(WLO.get_item_ij(k, d, (1,1), nu_sector)), # one band
								k, d)

		end 

	end # sector  

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function set_results_two!(results::AbstractDict,
													nk::Int, k0::Real,
													data)::Nothing

	for (dir1,d) in enumerate(data)
					
		set_results_one!(results, nk, k0, dir1, d) 

	end 

end  


function set_results_one!(results::AbstractDict, nk::Int, k0::Real,
													dir1::Int,
													(eigW1_occup,nus2pm)::Tuple{
																<:AbstractVector{<:AbstractArray},
																<:AbstractVector{<:AbstractArray{<:Real,3}}
																			},
#													(eigW1_unocc,eta2pm)::Tuple{
#																<:AbstractVector{<:AbstractArray},
#																<:AbstractVector{<:AbstractArray{<:Real,3}}
#																			},
													)::Nothing

	if haskey(results,"WannierGap")

		gap = WLO.WannierGap_fromSubspaces(eigW1_occup)
	
		setindex!(results["WannierGap"], gap, dir1)

	end  


	if haskey(results, "WannierBands1") 

		cp_wcc_to_results(results, "WannierBands1", 
											WLO.nuPlusMinus_fromSubspaces(eigW1_occup),
											dir1, nk)

	end 

	if haskey(results, "WannierBands2")

		cp_wcc_to_results(results, "WannierBands2", nus2pm, 3-dir1, nk)

	end 

	return 

end 


#===========================================================================#
#
#
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

	MBtime = braiding_time(P)
	nk = nr_kPoints(P)
	k0 = kPoint_start(P)
	symms = preserved_symmetries(P)
	strength = perturb_strength(P)
	parallel=nprocs()>=4


	results = init_results(nk, get_target(target; kwargs...))

#	all(isauxfile, keys(results)) && return results


	perturb1 = get_perturb_on_mesh(P,3268) # strength 1 so far, seed given


	psi = MB.get_psiH(MBtime, nk, k0, perturb1, strength)

	data = get_data(psi, results) 

	set_results_two!(results, nk, k0, data) 


	return results

end 











































































































#############################################################################
end # module CalcWLO

