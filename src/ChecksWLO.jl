module ChecksWLO
############################################################################# 


import LinearAlgebra 

import myLibs: ReadWrite, Utils 


import myLibs.Parameters:UODict 

import ..FILE_STORE_METHOD, ..WLO, ..MB


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


usedkeys()::Vector{Symbol} = [
						:braiding_time,
						:nr_kPoints,
						:kPoint_start,
						:preserved_symmetries,
						:nr_perturb_strength,
						:max_perturb_strength,
						:nr_perturb_instances,
						]

function usedkeys(P::UODict)::Vector{Symbol} 

	uk = usedkeys() 
	
	all_symms_preserved(P) && setdiff!(uk, [
								:nr_perturb_strength,
								:max_perturb_strength,
								:nr_perturb_instances
								]
					 )

	return uk

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function perturb_strengths(P::UODict)::Vector{Float64}

	LinRange(0,P[:max_perturb_strength],P[:nr_perturb_strength])

end  

function nr_perturb_instances(P::UODict)::Int 

	P[:nr_perturb_instances]

end  

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

function perturb_strengths_trials(P::UODict
																 )::Tuple{Vector{Float64}, 
																					Vector{Matrix{ComplexF64}} }
	
	all_symms_preserved(P) && return ([0.0], [zeros(ComplexF64,4,4)])
	
	randpert = rand(ComplexF64,4,4,nr_perturb_instances(P))

	return (perturb_strengths(P), collect.(eachslice(randpert,dims=3)))

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

braiding_theta(t::Real)::Float64  = t*2pi 

braiding_theta(P::UODict)::Float64 = braiding_theta(P[:braiding_time])



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function get_psiH(theta::Real, 
												 n::Int, k0::Real,
										 symms::AbstractString, 
										 args...
												 )

	if all_symms_preserved(symms)

			WLO.psiH_on_mesh(n, k0, MB.H,  MB.bsss_cycle(theta))
			
	else 

		WLO.psiH_on_mesh(n, k0, MB.perturbedH,  MB.bsss_cycle(theta), 
										 get_perturb(symms, args...))

	end 

end 
 
function psi_occup(psi::AbstractArray{ComplexF64,4}
									)::AbstractArray{ComplexF64,4}

	WLO.psi_sorted_energy(psi; halfspace=true, occupied=true)  

end 
function psi_unocc(psi::AbstractArray{ComplexF64,4}
									)::AbstractArray{ComplexF64,4}

	WLO.psi_sorted_energy(psi; halfspace=true, occupied=false)

end  

function get_data(psiH::AbstractArray{ComplexF64,4})

	(get_data(psiH, 1), get_data(psiH,2))

end 


function get_data(psiH::AbstractArray{ComplexF64,4}, dir1::Int,
									)

	(get_data(psiH, true, dir1), get_data(psiH, false, dir1))

end 

function get_data(psiH::AbstractArray{ComplexF64,4}, 
									occupied::Bool,
									dir1::Int,
									)

	psi = WLO.psi_sorted_energy(psiH; halfspace=true, occupied=occupied) 

	eigW1 = WLO.Wannier_subspaces_on_mesh_loseW(WLO.wlo1_on_mesh(dir1, psi))

	return (eigW1, WLO.wcc2mesh_fromSubspaces1(3-dir1, eigW1, psi))

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function init_results(L::Int)::Dict{String,Array}

	results = Dict{String,Array}()

	for eq in ("D123","D127a")

		results[eq] = zeros(L,2,2,2)

		results[eq*"_legend"] = ["dir1","sector","stat"]

	end 


	for eq in ("D110","D111","D125")

		results[eq] = zeros(L,2,2)
		results[eq*"_legend"] = ["dir1", "stat"]  

	end 

	for eq in ("WannierGap","D30","D48")
		
		results[eq] = zeros(L,2)

		results[eq*"_legend"] = ["stat"]

	end 

	results["WannierGap_legend"] = ["dir1"]

	return results

end 

function set_results_two!(results::AbstractDict,
													nk::Int, k0::Real,
													trial::Int,i_ps::Int,
													data)::Nothing

	for (dir1,d) in enumerate(data)
					
		set_results_one!(results, nk, k0, trial, i_ps, dir1, d...)

	end 

	return 


end 

function set_results_one!(results::AbstractDict, nk::Int, k0::Real,
													trial::Int,i_ps::Int,dir1::Int,
													(eigW1_occup,nus2pm)::Tuple{
																			<:AbstractVector{<:AbstractArray},
																			<:AbstractVector{<:AbstractArray}
																			},
													(eigW1_unocc,eta2pm)::Tuple{
																			<:AbstractVector{<:AbstractArray},
																			<:AbstractVector{<:AbstractArray}
																			},
													)

#function set_results_one!(results::AbstractDict,
#													trial::Int,i_ps::Int,dir1::Int,
#													eigW1_occup::AbstractVector{<:AbstractArray},
#													nus2pm::AbstractVector{<:AbstractArray},
#													eigW1_unocc::AbstractVector{<:AbstractArray},
#													eta2pm::AbstractVector{<:AbstractArray}
#													)
#

#	WannierGap 					
	gap = min(WLO.WannierGap_fromSubspaces(eigW1_occup),
						WLO.WannierGap_fromSubspaces(eigW1_unocc))
		
	WLO.run_m1!(results["WannierGap"], trial, gap, i_ps, dir1)

	
	p1occup = WLO.polariz_fromSubspaces!(eigW1_occup)
	p1unocc = WLO.polariz_fromSubspaces!(eigW1_unocc)

# eigW1[2] are overwritten

#	D110: p1occup-p1unocc ==0
	WLO.run_m1!(results["D110"],trial,
							WLO.wcc_stat_axpy(-1,p1occup,p1unocc),
							i_ps,dir1,:)
	
	

#	D111: p1occup+p1unocc ==0
	WLO.run_m1!(results["D111"],trial,
							WLO.wcc_stat_axpy(1,p1occup,p1unocc),
							i_ps,dir1,:)




	

	p1occup_ = sum(nus2pm)
	p1unocc_ = sum(eta2pm) 


#	D127a: polariz = sum over Wannier sectors  (occ=1/unocc=2)
	WLO.run_m1!(results["D127a"], trial, 
							WLO.wcc_stat_axpy!(-1,p1occup_,p1occup), 
						 i_ps,dir1,1,:) 

# p1occup has been overwritten

	WLO.run_m1!(results["D127a"], trial, 
							WLO.wcc_stat_axpy!(-1,p1unocc_,p1unocc), 
						 i_ps,dir1,2,:) 

#	p1unocc has been overwritten  

#	D125: total polarization zero 
	WLO.run_m1!(results["D125"], trial,
							WLO.wcc_stat_axpy!(1,p1occup_,p1unocc_),
						 i_ps,dir1,:)

#	p1unocc_ has been overwritten

	for sector in 1:2 # plus and minus 

#	D123: nu2pm == eta2pm 
		WLO.run_m1!(results["D123"], trial,
							 WLO.wcc_stat_axpy!(-1, nus2pm[sector], eta2pm[sector]), 
							 i_ps, dir1, sector, :) 

#	eta2pm has been overwritten

	end # sector 



	if dir1==1 
#		D30 for (d1,d2)==(1,2), one Wannier band per sector

		nu2_minus_mkx = view(nus2pm[2], 1,
												 WLO.ind_minusk.(axes(nus2pm[2],2),nk,k0),
												 :)

		WLO.run_m1!(results["D30"], trial,
							 WLO.wcc_stat_axpy!(-1,nus2pm[1],nu2_minus_mkx),
							 i_ps,:)

	elseif dir1==2
#	D48 for (d1,d2)==(2,1)
#		
		WLO.run_m1!(results["D48"], trial,
							 WLO.wcc_stat_axpy!(-1,nus2pm...),
							 i_ps,:)


	end #if 
	
 # nus2pm[2] has been overwritten

	return 

end 

function get_perturb(
										 symms::AbstractString, 
										 seed::AbstractMatrix{ComplexF64},
										 ps::Real,
										 )::Matrix{ComplexF64}

	all_symms_preserved(symms) && return zeros(ComplexF64, size(seed)...)

	LinearAlgebra.normalize!(MB.symmetrize(seed+seed', symms)) .*= ps 

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function Compute(P::UODict; get_fname::Function, target=nothing,
								 kwargs...)::Dict 

	theta = braiding_theta(P)
	nk = nr_kPoints(P)
	k0 = kPoint_start(P)
	symms = preserved_symmetries(P)


	strengths, trials = perturb_strengths_trials(P)

	results = init_results(length(strengths))


	for (j,perturb0) in enumerate(trials) # don't combine for loops! Order imp.

		for (i_ps,ps) in enumerate(strengths) 

			data = get_data(get_psiH(theta, nk, k0, symms, perturb0, ps))

			set_results_two!(results, nk, k0, j, i_ps, data)

		end 
	end 





	ReadWrite.Write_PhysObs(get_fname(P), FILE_STORE_METHOD, results)

	isnothing(target) && return results 

	return Utils.dict_keepkeys(results, 
														 vcat(([t,t*"_legend"] for t in vcat(target))...)
														 )

end  




#===========================================================================#
#
# adapted from 
# /media/tudor/Tudor/Work/2020_Snake-states/SnakeStates/Helpers/GF.jl
#
# Minimally edited: prep_obs -> vcat 
#
#---------------------------------------------------------------------------#


function FoundFiles(P; observables=nothing, target=nothing,
										get_fname::Function, kwargs...
									 )::Bool

	# kwargs for get_target

	FoundFiles0(get_fname(P),
							vcat(isnothing(observables) ? String[] : observables,
									 isnothing(target) ? String[] : target),
							)

end 





function FoundFiles0(Obs_fname::Function)::Function 

	function FFO(observables::Union{AbstractString, 
																	AbstractVector{<:AbstractString}})::Bool

		FoundFiles0(Obs_fname, observables)

	end 

end 


function FoundFiles0(Obs_fname::Function, observables::Union{AbstractString, 
																									AbstractVector{<:AbstractString}})::Bool

	ReadWrite.FoundFiles_PhysObs(Obs_fname, 
															 vcat(observables),
															 FILE_STORE_METHOD)


end





function Read(P::UODict; target=nothing,
							get_fname::Function, kwargs...
						 )::Dict 

	Read0(get_fname(P))

end
	



function Read0(Obs_fname::Function, 
												 ::Nothing=nothing,
												 ::Nothing=nothing)::Dict{String,Any}

	files_available = [split(fn,".")[1] for fn in cd(readdir, Obs_fname())]

	return ReadWrite.Read_PhysObs(Obs_fname, files_available, FILE_STORE_METHOD)

end	


function Read0(Obs_fname::Function, 
												 observables::Union{AbstractString, 
																						<:AbstractVector{<:AbstractString}},
												 target::Union{AbstractString,
																			 <:AbstractVector{<:AbstractString}},
												 )::Dict{String,Any}

	Read0(Obs_fname, intersect(vcat(observables), vcat(target)))

end  

function Read0(Obs_fname::Function, ::Nothing,
												 observables::Union{AbstractString, 
																						<:AbstractVector{<:AbstractString}},
												 )::Dict{String,Any}

	Read0(Obs_fname, observables)

end 


function Read0(Obs_fname::Function,
												 observables::Union{AbstractString, 
																						<:AbstractVector{<:AbstractString}},
												 ::Nothing=nothing,
												 )::Dict{String,Any}

	isempty(observables) && return Dict{String,Any}()

	files_available = [split(fn,".")[1] for fn in cd(readdir, Obs_fname())]

	filter!(f->any(t->occursin(t,f), vcat(observables)), files_available)

	return ReadWrite.Read_PhysObs(Obs_fname, files_available, FILE_STORE_METHOD)

end 





































































































#############################################################################
end  # module ChecksWLO

