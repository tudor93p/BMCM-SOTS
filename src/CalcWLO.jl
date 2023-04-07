module CalcWLO 
#############################################################################

using Distributed 
import Random, LinearAlgebra, Statistics, Combinatorics

using OrderedCollections: OrderedDict 

import myLibs: ReadWrite, Utils

import myLibs.Parameters: UODict 

import ..FILE_STORE_METHOD

import ..WLO, ..Helpers 

import ..MB; MODEL=MB
#import ..BBH; MODEL=BBH 

import ..WLO: nr_kPoints, kPoint_start 
import ..Helpers: Symmetries 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function all_symms_preserved end 
function perturb_strength end 
zero_perturb_strength = <(1e-12)∘perturb_strength   


Dependencies = [MODEL,WLO]


usedkeys()::Vector{Symbol} = [
						:preserved_symmetries, 

						:perturb_strength,
						]
									 
function usedkeys(P::UODict)::Vector{Symbol} 

	uk = usedkeys() 
	
	all_symms_preserved(P) && setdiff!(uk, [ :perturb_strength, ]) 

	zero_perturb_strength(P) && setdiff!(uk, [:preserved_symmetries, :perturb_strength])

#	(all_symms_preserved(P)||zero_perturb_strength(P))&& return Symbol[]

	return uk

end 

calc_observables = ["WannierBands1","WannierBands2"] # no "WannierGap"



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




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

	all_symms_preserved(P) ? get(P, :perturb_strength, 0) : P[:perturb_strength]

end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function fn_legend(obs::AbstractString)::String 

	islegend(obs) ? obs : obs*"_legend"

end  
function fn_nolegend(obs_leg::AbstractString)::String 

	split(obs_leg,"_legend")[1]

end  

function islegend(obs::AbstractString)::Bool 
	
	s = "_legend" 

	return length(obs)>length(s) && obs[end-length(s)+1:end]==s 

end 

#function fn_kWCC(obs::AbstractString)::String 
#
#	(occursin("WannierBands",obs) && !iskWCC(obs)) ? obs*"_k" : obs 
#	
#end  
#
#function fn_no_kWCC(obs_leg::AbstractString)::String 
#
#	split(obs_leg,"_k")[1]
#
#end  
#
#function iskWCC(obs::AbstractString)::Bool 
#
#	s = "_k"
#
#	return length(obs)>length(s) && obs[end-length(s)+1:end]==s 
#
#end 
#


function xxlabel()::Vector{String}

	["ks","klabels"]

end   

#
#isauxfile = isxxlabel = in(xxlabel())
#
function xxlabel(data::AbstractDict)::Tuple{<:AbstractVector{<:Real},
																						<:AbstractVector{<:AbstractString}
																						}

	Tuple(getindex(data,k) for k=xxlabel())

end 




prep_obs(::Nothing=nothing)::Vector{String} = String[]

function prep_obs(obs::AbstractString)::Vector{String} 

	obs1 = fn_nolegend(obs)

	obs1 in calc_observables || return prep_obs()

	return vcat(obs1, fn_legend(obs1), xxlabel()) 

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

function ylabels(obs::AbstractString, legend::AbstractString
								)::Vector{String}

	if legend=="dir"   

#		obs[end]=='1' && return ["\$\\nu_x\$","\$\\nu_y\$"]
		obs=="WannierBands1" && return ["\\overline{x}","\\overline{y}"]

#		obs[end]=='2' && return ["\\tilde{\\nu}_x","\\tilde{\\nu}_y"]
		obs=="WannierBands2" && return ["\\overline{x}_y","\\overline{y}_x"]

#		obs=="Polarization" && return ["p_x", "p_y"]
		
		obs=="WannierDensity" && return ["\\rho(R_y)", "\\rho(R_x)"]


	end 
	
	legend=="sector" && return ["+", "-"]

	error("Wrong input")

end 


function add_legend!(results::AbstractDict, obs::AbstractString,
										 legend::AbstractVector{<:AbstractString}
										 )::Nothing



	setindex!(results, 
						OrderedDict((leg=>ylabels(obs,leg) for leg=legend)...), 
						fn_legend(obs))

	return 

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function init_results(n::Int, obs::AbstractVector{<:AbstractString}
										 )::Dict{String,Any}

	if "WannierBands2" in obs 

		union!(obs, calc_observables)

#	elseif "WannierBands1" in obs 

#		union!(obs, ["WannierGap"])

	end  

	results = Dict{String,Any}("klabels" => ["k_y","k_x"])

	for w in ("WannierBands1","WannierBands2")

		in(w,obs) || continue 

		results[w] = zeros(n,2,2)

		add_legend!(results, w, ["dir","sector"])

	end 
#	for w in ("WannierGap",)
#		
#		in(w,obs) || continue 
#
#		results[w] = zeros(2) 
#
#	end 



	return results

end 

function init_results(n::Int, k0::Real,
											obs::AbstractVector{<:AbstractString}
													 )::AbstractDict{String,Any}

	results = init_results(n, obs) 

	results["ks"] = WLO.get_kij(n, k0; restricted=false)(1:n)
	
	@assert Set(keys(results))==Set(xxlabel())

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


#function symmetrize_on_mesh(seed::AbstractArray{ComplexF64,4},
#														args...; kwargs...
#														 )::Array{ComplexF64,4}
#
#	A = copy(seed)
#
#	symmetrize_on_mesh!(A, args...; kwargs...)
#
#	return A
#
#end  
#
#


function symmetrize_and_normalize!(pert::AbstractArray{ComplexF64,4}, 
																	 args...)::AbstractArray{ComplexF64,4}

	WLO.store_on_mesh!!(Symmetries.symmetrize_HC!, pert)

	MODEL.symmetrize_on_mesh!(pert, args...)

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
											) 

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


function cp_wcc_to_results(results::AbstractDict{<:AbstractString,<:Any},
													obs::AbstractString,
													wccs::Union{Tuple{<:T,<:T}, AbstractVector{<:T}
																		 } where T<:AbstractArray{<:Real,3},
													d::Int,
													n::Int,
													k0::Real
													) 

	for (sector,nu_sector) in enumerate(wccs)

		WLO.check_nu_k_dep(selectdim(nu_sector,1,1), d) # one band 

		for i=1:n 

			setindex!(results[obs],
								only(WLO.get_item_ij(i, 3-d, (1,1), nu_sector)), # one band
								i, d, sector)

		end 

	end # sector  

end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function set_results!(results::AbstractDict,
													nk::Int, k0::Real,
													data)::Nothing

	for (dir1,d) in enumerate(data)
					
		set_results_onedir!(results, nk, k0, dir1, d) 

	end 

end  


function set_results_onedir!(results::AbstractDict, nk::Int, k0::Real,
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

#	if haskey(results,"WannierGap")
#
#		gap = WLO.WannierGap_fromSubspaces(eigW1_occup)
#	
#		setindex!(results["WannierGap"], gap, dir1)
#
#	end  


	if haskey(results, "WannierBands1") 


		cp_wcc_to_results(results, "WannierBands1", 
											WLO.nuPlusMinus_fromSubspaces(eigW1_occup),
											dir1, nk, k0)

	end 

	if haskey(results, "WannierBands2")


		cp_wcc_to_results(results, "WannierBands2", nus2pm, 3-dir1, nk, k0)

	end 

	return 

end 


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
	#parallel=nprocs()>=2 # parallel 2 cores not worth it 

	results = init_results(nk, k0, get_target(target; kwargs...))

	perturb1 = get_perturb_on_mesh(P)

#
#	if !in(symms ,["None","All"])
#
#	symms2 = join.(Combinatorics.powerset(["P", "Mx", "Ct", "Tt"],1,2),"+") 
#
#	println()
#
#	@show symms 
#
#	println("Ham.: ",
#					all(MODEL.has_symm_on_mesh(MODELparams, symms, nk,k0)))
#
#	for s2 in symms2 
#
#		println("$s2:\t\t",
#						all(MODEL.has_symm_on_mesh(MODELparams, s2, nk,k0)))
#
#	end 
#
#
#	println("\nPert.: ",
#					all(MODEL.has_symm_on_mesh(perturb1, symms, nk,k0)))
#
#	for s2 in symms2 
#
#		println("$s2:\t\t",
#						all(MODEL.has_symm_on_mesh(perturb1, s2, nk,k0)))
#
#	end 
#
#	println()
#
#
#end 


################## test 
#
#	perturb2 = get_perturb_on_mesh("Ct", nk, k0, 1993)
#
##	perturb1 = perturb1*strength + perturb2*1e-10 
#
#	LinearAlgebra.axpby!(1e-10, perturb2, strength, perturb1)
#
#	strength = 1
#
#	@show LinearAlgebra.norm(perturb1)
#
#########################

	psi = MODEL.get_psiH(P, nk, k0, perturb1, strength)


	set_results!(results, nk, k0, get_data(psi, results))

	return results

end 



#===========================================================================#
#
# FoundFiles &  Read 
#
#---------------------------------------------------------------------------#


function FoundFiles(P::UODict; 
										target=nothing, get_fname::Function, kwargs...)::Bool

	FoundFiles0(get_fname(P), get_target(target; kwargs...))

end 



function FoundFiles0(Obs_fname::Function, 
										 observables::AbstractVector{<:AbstractString})::Bool

	!isempty(observables) && ReadWrite.FoundFiles_PhysObs(Obs_fname, observables, FILE_STORE_METHOD)

end





function Read(P::UODict; target=nothing, get_fname::Function, kwargs...)::Dict 

	Read0(get_fname(P), get_target(target; kwargs...))

end
	



function Read0(Obs_fname::Function, 
							 observables::AbstractVector{<:AbstractString}
												 )::Dict{String,Any}

	if isempty(observables)
		
		obs_found = [split(fn,".")[1] for fn in cd(readdir, Obs_fname())]

		return ReadWrite.Read_PhysObs(Obs_fname, obs_found, FILE_STORE_METHOD)

	else 

		return ReadWrite.Read_PhysObs(Obs_fname, observables, FILE_STORE_METHOD)

	end 

end 










































































































#############################################################################
end # module CalcWLO

