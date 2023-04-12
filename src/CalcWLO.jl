module CalcWLO 
#############################################################################

#using Distributed  

import Random, LinearAlgebra, Statistics, Combinatorics, Optim

using OrderedCollections: OrderedDict 

import myLibs: ReadWrite, Utils, SignalProcessing 

import myLibs.Parameters: UODict 

import ..FILE_STORE_METHOD

import ..WLO, ..Helpers 

import ..MB; MODEL=MB
#import ..BBH; MODEL=BBH 

import ..WLO: nr_kPoints, kPoint_start 
import ..Helpers: Symmetries, AdaptiveMesh 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function all_symms_preserved end 
function perturb_strength end 
function adaptive_kMesh end 
zero_perturb_strength = <(1e-12)∘perturb_strength   





Dependencies = [MODEL,WLO]


usedkeys()::Vector{Symbol} = [

						:kMesh_type, :kMesh_model,

						:preserved_symmetries, 

						:perturb_strength,
						]
									 
function usedkeys(P::UODict)::Vector{Symbol} 

	uk = usedkeys() 
	
	adaptive_kMesh(P) || setdiff!(uk, [:kMesh_type, :kMesh_model])

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



function adaptive_kMesh(P::UODict)::Bool 

	haskey(P, :kMesh_type) && lowercase(P[:kMesh_type]) == "adaptive"

end 


function kMesh_model(P::UODict)::Function 

	@assert adaptive_kMesh(P)

	n = P[:kMesh_model]

	getf = getproperty(AdaptiveMesh, Symbol("f_"*n))

	getp = getproperty(AdaptiveMesh, Symbol(n*"_par_from_vals"))

	return getf∘getp 

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
function xxlabel(data::AbstractDict)::Tuple{<:AbstractMatrix{<:Real},
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

	results["ks"] = repeat(WLO.get_kij(n, k0; restricted=false)(1:n),
												 inner=(1,2))

	@assert issubset(xxlabel(), keys(results))

	return results

end  
function init_results(n::Int,
											(kx,ky)::AbstractVector{<:AbstractVector{<:Real}},
											obs::AbstractVector{<:AbstractString}
										 )::Dict{String,Any}

	results = init_results(n, obs) 

	results["ks"] = hcat(ky,kx) # consistent: klabels=["k_y","k_x"],W=zeros(n,2)

	@assert issubset(xxlabel(), keys(results))

	return results

end   



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function sum_kSteps(ks::AbstractVector{<:Real},
										nk::Int
										 )::Float64

	check_vec_len(nk,ks)

	@assert length(ks)==nk 

	return abs(ks[end]-ks[1])  

end 

function sum_kSteps_dist2pi!(
											 gaps::AbstractVector{Float64},
											 ks::AbstractVector{Float64},
												 nk::Int, 
												 args...
													 )::Function 

	function ssi2p(alpha::Union{AbstractVector{<:Real},<:Real})::Float64 
	
		fill_gaps_ks!(gaps, ks, nk, args..., only(alpha))
	
#		println(only(alpha),"\t",sum_kSteps(ks))

		return abs(2pi-sum_kSteps(ks,nk))

	end 

end  




function init_gaps_ks(nk::Int)::NTuple{2,Vector{Float64}}

	Tuple(WLO.init_storage(Float64, nk+1) for i=1:2)

end 

function check_vec_len(nk::Int, vectors::AbstractVector{<:Real}...)

	for v in vectors 
		@assert WLO.nr_kPoints_from_mesh1(v)==nk+1
	end 

end 




function fill_gaps_ks(nk::Int, 
											args...
											)::NTuple{2,Vector{Float64}}

	fill_gaps_ks!(init_gaps_ks(nk)..., nk, args...)


end 


function next_kPoint(get_dk_for_gap::Function,
										 gap::Real,
										 k_prev::Real,
										br_args...)

	dk = get_dk_for_gap(gap)

	@assert dk>0 "Decrease minimum expected gap"

	WLO.next_kPoint(k_prev, 
									AdaptiveMesh.bound_rescale_kStep(dk, br_args...)
									)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function wrapper_get_gap(get_gap_at_k::Union{Function,
																						 SignalProcessing.Dierckx.Spline1D},
												 k::Real
												 )::Float64


	get_gap_at_k(k)

end 
function wrapper_get_gap(
												 (get_gap_at_k,data_gap)::Tuple{Function,<:Any},
												 k::Real
												 )::Float64

	get_gap_at_k(data_gap, k)

end 


function fill_gaps_ks!(
											 gaps_::AbstractVector{Float64},
											 ks_::AbstractVector{Float64},

												 nk::Int, k0::Real,
																	unif_kij::Function,
																	uniq_kinds::AbstractVector{Int},
																	ind_minusk::Function,
																	
																	get_gap_at_k::Union{<:Tuple{Function,<:Tuple}, Function, SignalProcessing.Dierckx.Spline1D},

											get_dk_for_gap::Function,

											br_args...

											)::Tuple{<:AbstractVector{Float64},
															 <:AbstractVector{Float64}
															 }

	check_vec_len(nk, gaps_, ks_)

	gaps = view(gaps_, uniq_kinds)
	ks = view(ks_, uniq_kinds)

	setindex!(ks, unif_kij(uniq_kinds[1])[1], 1)

	for i=2:length(uniq_kinds)

		setindex!(gaps, wrapper_get_gap(get_gap_at_k, ks[i-1]), i-1)

		setindex!(ks, 
							next_kPoint(get_dk_for_gap, gaps[i-1], ks[i-1], br_args...),
							i) 

	end  
		
	setindex!(gaps, wrapper_get_gap(get_gap_at_k, ks[end]), length(gaps)) 

	for i in uniq_kinds 

		j = ind_minusk(i)

		setindex!(ks_, -ks_[i], j)

		setindex!(gaps_, gaps_[i], j)

	end 
	
	
	#	ks_[nk] and gaps_[nk] unaccessed so far 
		
	setindex!(ks_, 
						next_kPoint(get_dk_for_gap, gaps_[nk-1], ks_[nk-1], br_args...), 
						nk)

	setindex!(gaps_, wrapper_get_gap(get_gap_at_k, ks_[nk]), nk)



	return gaps, ks

end 




function find_rescaling_factor_dk(nk::Int, args...; kwargs...)

	find_rescaling_factor_dk!(init_gaps_ks(nk)..., nk, args...; kwargs...)

end 


function find_rescaling_factor_dk!(
											 gaps::AbstractVector{Float64},
												ks::AbstractVector{Float64}, 

												 nk::Int, k0::Real, 

												 get_gap_at_k,

											get_dk_for_gap::Function,

																	bounds::AbstractVector{<:Real};
																	optim_tol::Float64=1e-8

								

														 )#::Tuple{Float64,Vector{Float64}}

#@assert length(gaps)==length(ks)==nk  
	check_vec_len(nk, gaps, ks)

	bounds_new = AdaptiveMesh.verify_dk_bounds(bounds, nk) 

	uniq_kinds = WLO.uniqueInds_kMirror(nk,k0)
	unif_kij = WLO.get_kij(nk,k0) #restricted=false)
	ind_minusk = WLO.ind_minusk(nk,k0)


	fill_gaps_ks!(gaps, ks, nk, k0, 
																	unif_kij,
																	uniq_kinds,
																	ind_minusk,

								get_gap_at_k, get_dk_for_gap, 
								bounds_new,
																	)


	sksd2p = sum_kSteps_dist2pi!(
								gaps, ks, 
								nk, k0,
																	unif_kij,
																	uniq_kinds,
																	ind_minusk,


								get_gap_at_k,
								get_dk_for_gap, 
								bounds_new,
								)
	sol = Optim.optimize(sksd2p, [1.0],
											 Optim.Options(g_tol = optim_tol,
																		 f_abstol=optim_tol,
																		 ))

	alpha = only(Optim.minimizer(sol))

	get_dk_for_gap_new = AdaptiveMesh.bound_rescale_kStep(get_dk_for_gap, bounds_new, alpha)

	fill_gaps_ks!(gaps, ks, nk, k0, 
																	unif_kij,
																	uniq_kinds,
																	ind_minusk,
								get_gap_at_k, get_dk_for_gap_new)



#	correct_sum_kSteps!(ks, nk, uniq_kinds, ind_minusk)

	return (gaps, ks), get_dk_for_gap_new 

end 





























#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

"""
sector = (occupied,dir1)
nks = (rough_nk, decent_nk, dense_nk)

Compute the gap of W_dir1 for k_dir1 uniform and k_dir2 adaptive.

Find the parameter vector ks_dir2 which is more dense for small gaps
"""
function threestep_find_ks(Hdata, (nks,k0), 
													 sector,
						(extrema_gaps, extrema_dk),
						model::Function
						)

	@assert length(nks)==3 && issorted(nks)


	# --- rough estimation of gap(k) in a halfspace --- #
	
	data_1 = WLO.pack_data_gap(Hdata, (nks[1],k0), sector)

	gaps_1, ks_1 = [WLO.calc_gap!(data_1,k) for k=WLO.uniqueKs_kMirror(nks[1],k0)]



	min_gap = min(extrema_gaps[1], minimum(gaps_1))
	
	dk_from_gap_1 = model([min_gap, extrema_gaps[2]], extrema_dk)
	



	# --- linearly scale 'model' such that ks span a halfspace --- #

	gap_at_k_2 = (WLO.calc_gap!,WLO.pack_data_gap(Hdata, (nks[2],k0), sector)) 

	(gaps_2, ks_2),dk_from_gap_2 = find_rescaling_factor_dk(
																													nks[2], k0, 
																													gap_at_k_2,
																							dk_from_gap_1, extrema_dk;
																							optim_tol=1e-8
																							)


	# --- scale again, full precision --- # 

	if ks_2[2]<ks_2[1] 

		reverse!(ks_2)
		reverse!(gaps_2)

	end 

	gap_at_k_3 = SignalProcessing.Interp1D(ks_2, gaps_2, 3) 
	
	(gaps_3, ks_3),dk_from_gap_3 = find_rescaling_factor_dk(nks[3], k0, 
																													gap_at_k_3,
																							dk_from_gap_2, extrema_dk;
																							optim_tol=1e-14
																						)

	return (gaps_3, ks_3, dk_from_gap_3,
					(nks[3],k0),
					sector)

end 	


function kxy_from_threestep(
														 (gaps_1, precomp_ks_1, get_dk_1,
															(nk_1,k0_1), (occ_1, perp_dir_1)
															),
														 (gaps_2, precomp_ks_2, get_dk_2,
															(nk_2,k0_2), (occ_2, perp_dir_2)
															),
														 ;
														 kwargs...
														 )::Vector{Vector{Float64}}


	@assert perp_dir_1!=perp_dir_2 

	perp_dir_1==2 ? [precomp_ks_1,precomp_ks_2] : [precomp_ks_2,precomp_ks_1]

end 



function calc_kxy_adaptive(Hdata, 
													 nk::Int,
													 k0::Real, 
													 step_vs_gap_model::Function 
													 )::Vector{Vector{Float64}}

	nks = [min(17,nk),min(71,nk),nk] 

	extrema_gaps_dk = ([0.02, 0.5],[1e-7, pi/15]) 

	out_1,out_2 = map(1:2) do dir 
		
		threestep_find_ks(Hdata, (nks,k0), (true,dir), extrema_gaps_dk, step_vs_gap_model) 

	end 

	return kxy_from_threestep(out_1, out_2)

end 

function calc_kxy_adaptive(P::UODict 
													 )::Vector{Vector{Float64}}

	@assert adaptive_kMesh(P) 

	# no perturbation to determine mesh   
	
	calc_kxy_adaptive(MODEL.get_args_psi(P),
										nr_kPoints(P),
										kPoint_start(P),
										kMesh_model(P),
													 )

end 




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
#									parallel::Bool=false
									)::Vector

#	(parallel ? pmap : map)
	map(Base.splat(WLO.get_wlo_data_mesh), get_data_args(psiH, results))

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
													nk::Int, #k0::Real,
													data)::Dict 

	for (dir1,d) in enumerate(data)
					
		set_results_onedir!(results, nk,  dir1, d) 

	end 
	
	return results 

end  


function set_results_onedir!(results::AbstractDict, nk::Int, #k0::Real,
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
											dir1, nk)

	end 

	if haskey(results, "WannierBands2")


		cp_wcc_to_results(results, "WannierBands2", nus2pm, 3-dir1, nk)

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

	Compute_(P, target, 
					 adaptive_kMesh(P) ? calc_kxy_adaptive(P) : kPoint_start(P);
					 kwargs...)

end  



function Compute_( P::UODict,
									 target,
											k::Union{<:Real,
															 <:AbstractVector{<:AbstractVector{<:Real}}
															 };
									kwargs...)::Dict{String,Any}


	nk = nr_kPoints(P)   

	results = init_results(nk, k, get_target(target; kwargs...))  

	psi = MODEL.get_psiH(P, nk, k, 
											 get_perturb_on_mesh(P), 
											 perturb_strength(P),
											 ) 

	return set_results!(results, nk, get_data(psi, results))

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

