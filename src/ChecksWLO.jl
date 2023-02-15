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

	preserved_symmetries(P)=="All"

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

function get_partial_data(theta::Real, perturb::AbstractMatrix{<:Number},
												 n::Int, k0::Real,
												 )
#	@show LinearAlgebra.norm(perturb)

	data = MB.get_wcc1_data(n, k0, MB.perturbedH,  
															MB.bsss_cycle(theta),
															perturb)

	data_gap = [MB.get_wcc1gap(data,d1) for d1=1:2]

	data_wcc2 = [WLO.wcc_stat!(MB.get_wcc2(data, 3-d1), MB.quantized_wcc2_values; dim=2)[1,:] for d1=1:2]

	return data_gap,data_wcc2

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function Compute(P::UODict; get_fname::Function, target=nothing,
								 kwargs...)::Dict 

	results = Dict{String,Any}("D110"=>[1,2,3])

########## 

	nk = nr_kPoints(P)
	k0 = kPoint_start(P)
	theta = braiding_theta(P)

	strengths = perturb_strengths(P) 

#	if all_symms_preserved(P)
#
#		println("* idling *")
#
#
#	else 


	nr_trials = nr_perturb_instances(P) 



	symms = preserved_symmetries(P)


	randpert = [p+p' for p in eachslice(rand(ComplexF64,4,4,nr_trials),dims=3)] 


	Y2Stat = zeros(length(strengths), 3, 2)

	selectdim(Y2Stat,2,3) .= 1 



	for (j,perturb0) in enumerate(randpert) 

		perturb = LinearAlgebra.normalize!(MB.symmetrize(perturb0, symms))
	
#		perturb *= maximum(strengths) 



#		data = (nprocs()<=3 ? map : pmap)(strengths) do ps 
#
#			get_partial_data(theta, ps*perturb, nk, k0)
#
#		end 

#
#		for (i_ps,(ps,(data_gap,data_wcc2))) in enumerate(zip(strengths,data))
		for (i_ps,ps) in enumerate(strengths)
#				data_gap,data_wcc2 = data[i_ps]

				data_gap,data_wcc2 = get_partial_data(theta, ps*perturb, nk, k0)

#				data_= WLO.MB.get_wcc1_data(WLO.MB.perturbedH,
#																		WLO.MB.bsss_cycle(theta),
#																		ps*perturb)
#	

##

#				fig1, Ax1,  = init_fig(theta; num=i_ps)
#
#
				for d1=1:2 
					
					WLO.run_m1!(Y2Stat, j, data_gap[d1], i_ps, 3, d1)
					
					for (i,w) in enumerate(data_wcc2[d1])
						
						WLO.run_m1!(Y2Stat, j, abs(w), i_ps, i, d1)

					end  

				end # d1 


		end # ps 
	end # trials 

	@show Y2Stat

	@show get_fname(P)("D110")

















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

