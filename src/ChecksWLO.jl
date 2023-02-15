module ChecksWLO
############################################################################# 

import myLibs: ReadWrite 

import myLibs.Parameters:UODict 

import ..FILE_STORE_METHOD, ..WLO 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


usedkeys = [:nr_kPoints,
						:preserved_symmetries,
						:nr_perturb_strength,
						:max_perturb_strength,
						:nr_perturb_instances,
						:kPoint_start,
						]


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





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function get_partial_data(theta::Real, perturb::AbstractMatrix{<:Number})
#	@show LinearAlgebra.norm(perturb)

	data = WLO.MB.get_wcc1_data(WLO.MB.perturbedH,
															WLO.MB.bsss_cycle(theta),
															perturb)

	data_gap = [WLO.MB.get_wcc1gap(data,d1) for d1=1:2]

	data_wcc2 = [WLO.wcc_stat!(WLO.MB.get_wcc2(data, 3-d1), WLO.MB.quantized_wcc2_values; dim=2)[1,:] for d1=1:2]

	return data_gap,data_wcc2

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function Compute(P::UODict; get_fname::Function, target=nothing,
								 kwargs...)::Dict 

	results = Dict{String,Any}("Hello"=>[1,2,3])

########## 

	perturb_strength = perturb_strengths(P) 

	nr_trials = nr_perturb_instances(P) 

	nk = nr_kPoints(P)


	randpert = [p+p' for p in eachslice(rand(ComplexF64,4,4,nr_trials),dims=3)] 


	Y2Stat = zeros(length(perturb_strength), 3, 2)

	selectdim(Y2Stat,2,3) .= 1 



















##########

	ReadWrite.Write_PhysObs(get_fname(P), FILE_STORE_METHOD, 
													results,
													# keys(results),
													)


	isnothing(target) && return results 

	return Utils.dict_keepkeys(results, vcat(target)) 

end  




#===========================================================================#
#
# adapted from 
# /media/tudor/Tudor/Work/2020_Snake-states/SnakeStates/Helpers/GF.jl
#
# Minimally edited: prep_obs -> vcat 
#
#---------------------------------------------------------------------------#


function FoundFiles(P; target=nothing,
										get_fname::Function, kwargs...
									 )::Bool

	# kwargs for get_target
	
	FoundFiles0(get_fname(P))
	

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

