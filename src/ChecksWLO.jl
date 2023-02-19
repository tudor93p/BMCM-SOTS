module ChecksWLO
############################################################################# 


import LinearAlgebra

using OrderedCollections: OrderedDict 

import myLibs: ReadWrite, Utils 


import myLibs.Parameters:UODict 

import ..FILE_STORE_METHOD, ..WLO, ..MB, ..Helpers



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

function get_data(psiH::AbstractArray{ComplexF64,4}, 
									results::AbstractDict)

	(get_data(psiH, 1, results), get_data(psiH,2, results))

end 


function get_data(psiH::AbstractArray{ComplexF64,4}, dir1::Int,
									results::AbstractDict
									)

	(get_data(psiH, true, dir1, results), get_data(psiH, false, dir1, results))

end 

obs_batch_1 = ["D110", "D111", "WannierGap"]
obs_batch_2 = ["D30", "D48","D123", "D125", "D127a"]  

function get_data(psiH::AbstractArray{ComplexF64,4}, 
									occupied::Bool,
									dir1::Int,
									results::AbstractDict
									)




	psi = WLO.psi_sorted_energy(psiH; halfspace=true, occupied=occupied) 

	eigW1 = WLO.Wannier_subspaces_on_mesh_loseW(WLO.wlo1_on_mesh(dir1, psi))

	if any(in(keys(results)), obs_batch_2)

		return (eigW1, WLO.wcc2mesh_fromSubspaces1(3-dir1, eigW1, psi))

	else 

		return (eigW1, fill(zeros(0,0,0),0))

	end

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

	occursin("_legend",obs)

end 


function add_legend!(results::AbstractDict, obs::AbstractString,
										 legend::AbstractVector{<:AbstractString}
										 )::Nothing


	setindex!(results, 
						OrderedDict((leg=>ylabels(obs,leg) for leg=legend)...), 
						fn_legend(obs))

	return 

end 


function get_label(c::CartesianIndex, args...)::String

	get_label(c.I, args...)

end 

function get_label(I::Tuple{Vararg{Int}}, legend::OrderedDict)::String

	get_label(I, values(legend))

end 


function get_label(I::Tuple{Vararg{Int}}, legend::Union{<:AbstractVector{<:AbstractVector},
																												<:Base.ValueIterator})::String

	#	[string(leg,": \$",ylab[i],"\$") for (i,(leg,ylab)) in zip(c.I,pairs(legend))]

	isempty(I) && return ""

	return "\$"*join(map(Base.splat(getindex), zip(legend, I)),",")*"\$"

end 



function xxlabel()::Vector{String}

	["x","xlabel"]

end  

function xxlabel(data::AbstractDict)::Tuple{<:AbstractVector,
																						<:AbstractString}

	Tuple(getindex(data,k) for k=xxlabel())

end 


isxxlabel = in(xxlabel())

function isauxfile(obs::AbstractString)::Bool 

	for f in (isxxlabel, islegend)

		f(obs) && return true 

	end 

	return false 

end 


function obs_CI(Y::AbstractArray{<:Real,N}, 
								obs_i::Int, ax::AbstractVector{Int}=2:N 
							 )::CartesianIndex where N
		
	CartesianIndices(axes(Y)[ax])[min(end,obs_i)]

end 


function get_data_and_lab(Y::AbstractArray, 
													legend::OrderedDict,
													obs_i::Int)::Tuple{<:AbstractVector,
																						 <:AbstractString}

	c = obs_CI(Y, obs_i)

	return view(Y, :, c), get_label(c, legend)

end  


function get_data_and_lab(data::AbstractDict,
													obs::AbstractString,
													obs_i::Int)::Tuple{<:AbstractVector,
																						 <:AbstractString}

	get_data_and_lab(data[fn_nolegend(obs)], data[fn_legend(obs)], obs_i)

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function init_results(strengths::AbstractVector{<:Real},
											obs::AbstractVector{<:AbstractString}
										 )::Dict{String,Any}


	if any(in(obs), obs_batch_2)

		union!(obs, obs_batch_1, obs_batch_2)

	elseif any(in(obs), obs_batch_1)

		union!(obs, obs_batch_1)

	end 

	results = Dict{String,Any}("x" => strengths,
														 "xlabel"=>"Perturbation strength"
															)

	L = length(strengths)

	for eq in ("D123","D127a",)

		in(eq,obs) || continue 

		results[eq] = zeros(L,2,2,2)

		add_legend!(results, eq, ["dir1","sector","stat"])

	end 


	for eq in ("D110","D111","D125",) 

		in(eq,obs) || continue 

		results[eq] = zeros(L,2,2)

		add_legend!(results, eq, ["dir1", "stat"])

	end 

	for eq in ("D30","D48",) 

		in(eq,obs) || continue 
		
		results[eq] = zeros(L,2)

		add_legend!(results, eq, ["stat"])

	end 

	for eq in ("WannierGap",)

		in(eq,obs) || continue 
		
		results[eq] = zeros(L,2)

		add_legend!(results, eq, ["dir1"])

	end 

	return results

end 



function ylabels(obs::AbstractString, legend::AbstractString
								)::Vector{String}

	legend=="dir1" && return ["x", "y"]
	
	legend=="stat" && return ["\\mu","\\sigma"]

	if legend=="sector" 

		obs=="D123" && return ["+", "-"]

#		obs=="D127a" && return ["\\\\text{occup}","\\\\text{unocc}"]
		obs=="D127a" && return ["\\mathrm{occup}","\\mathrm{unocc}"]
#
# tex command \text{...} not understood in plots 
#

	end  

	error("Wrong input")

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



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
													)::Nothing

#	WannierGap 					
#
	if haskey(results,"WannierGap")

		gap = min(WLO.WannierGap_fromSubspaces(eigW1_occup),
							WLO.WannierGap_fromSubspaces(eigW1_unocc))
	
	
		WLO.run_m1!(results["WannierGap"], trial, gap, i_ps, dir1)

	end 



	# -------------------------------------------- #

	if any((haskey(results,k) for k in ("D110","D111","D127a","D125")))


		p1occup = WLO.polariz_fromSubspaces!(eigW1_occup)
		p1unocc = WLO.polariz_fromSubspaces!(eigW1_unocc)
	
	# eigW1[2] are overwritten

		if haskey(results,"D110")
	
			#	D110: p1occup-p1unocc ==0
				WLO.run_m1!(results["D110"],trial,
										WLO.wcc_stat_axpy(-1,p1occup,p1unocc),
										i_ps,dir1,:)
			
		end 
		
		if haskey(results,"D111")

		#	D111: p1occup+p1unocc ==0
			WLO.run_m1!(results["D111"],trial,
									WLO.wcc_stat_axpy(1,p1occup,p1unocc),
									i_ps,dir1,:)
		
		end 	
	
	
		if haskey(results,"D127a")|haskey(results,"D125")		
	
			p1occup_ = sum(nus2pm)
			p1unocc_ = sum(eta2pm) 
		
			if haskey(results,"D127a")
	
			#	D127a: polariz = sum over Wannier sectors  (occ=1/unocc=2)
				WLO.run_m1!(results["D127a"], trial, 
										WLO.wcc_stat_axpy!(-1,p1occup_,p1occup), 
									 i_ps,dir1,1,:) 
			
			# p1occup has been overwritten
			
				WLO.run_m1!(results["D127a"], trial, 
										WLO.wcc_stat_axpy!(-1,p1unocc_,p1unocc), 
									 i_ps,dir1,2,:) 
			
			#	p1unocc has been overwritten  
	
			end 
	
			if haskey(results,"D125")
	
			#	D125: total polarization zero 
				WLO.run_m1!(results["D125"], trial,
										WLO.wcc_stat_axpy!(1,p1occup_,p1unocc_),
									 i_ps,dir1,:)
			
			#	p1unocc_ has been overwritten
			end 	

		end 

	end 

	# -------------------------------------------- #


	# -------------------------------------------- #
	for sector in 1:2 # plus and minus 

#	D123: nu2pm == eta2pm 
		haskey(results,"D123") && WLO.run_m1!(results["D123"], trial,
							 WLO.wcc_stat_axpy!(-1, nus2pm[sector], eta2pm[sector]), 
							 i_ps, dir1, sector, :) 

#	eta2pm has been overwritten

	end # sector 



	if dir1==1 && haskey(results,"D30")
#		D30 for (d1,d2)==(1,2), one Wannier band per sector

		nu2_minus_mkx = view(nus2pm[2], 1,
												 WLO.ind_minusk.(axes(nus2pm[2],2),nk,k0),
												 :)

		WLO.run_m1!(results["D30"], trial,
							 WLO.wcc_stat_axpy!(-1,nus2pm[1],nu2_minus_mkx),
							 i_ps,:)

	elseif dir1==2 && haskey(results,"D48")
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

	theta = braiding_theta(P)
	nk = nr_kPoints(P)
	k0 = kPoint_start(P)
	symms = preserved_symmetries(P)


	strengths, trials = perturb_strengths_trials(P)


	results = init_results(strengths, get_target(target; kwargs...))

#	return Dict{String,Any}(k=>v for (k,v)=pairs(results) if isauxfile(k)) 


	all(isauxfile, keys(results)) && return results

	for (j,perturb0) in enumerate(trials) # don't combine for loops! Order imp.

		for (i_ps,ps) in enumerate(strengths) 

			data = get_data(get_psiH(theta, nk, k0, symms, perturb0, ps), results)

			set_results_two!(results, nk, k0, j, i_ps, data)

		end 
	end 

	return results

end 





#===========================================================================#
#
# adapted from 
# /media/tudor/Tudor/Work/2020_Snake-states/SnakeStates/Helpers/GF.jl
#
# Minimally edited: prep_obs -> vcat 
#
#---------------------------------------------------------------------------#


prep_obs(::Nothing=nothing)::Vector{String} = String[]

function prep_obs(obs::AbstractString)::Vector{String} 
	
	vcat(fn_nolegend(obs),fn_legend(obs),xxlabel())

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





function FoundFiles(P; target=nothing, get_fname::Function, kwargs...)::Bool

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


#function Read0(Obs_fname::Function,
#												 observables::Union{AbstractString, 
#																						<:AbstractVector{<:AbstractString}},
#												 ::Nothing=nothing,
#												 )::Dict{String,Any}
#
#	isempty(observables) && return Dict{String,Any}()
#
#	files_available = [split(fn,".")[1] for fn in cd(readdir, Obs_fname())]
#
#	filter!(f->any(t->occursin(t,f), vcat(observables)), files_available)
#
#	return ReadWrite.Read_PhysObs(Obs_fname, files_available, FILE_STORE_METHOD)
#
#end 





































































































#############################################################################
end  # module ChecksWLO

