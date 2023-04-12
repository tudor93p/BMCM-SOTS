module ChecksWLO
############################################################################# 


using Distributed 

import LinearAlgebra, Combinatorics, Random 

using OrderedCollections: OrderedDict 

import myLibs: ReadWrite, Utils, SignalProcessing


import myLibs.Parameters:UODict 

import ..FILE_STORE_METHOD, ..WLO, ..Helpers


import ..Helpers: Symmetries

import myPlots

import ..WLO: nr_kPoints, kPoint_start  

import ..CalcWLO
import ..CalcWLO: preserved_symmetries, all_symms_preserved, MODEL

#get_perturb_on_mesh_, 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

Dependencies = [CalcWLO]


usedkeys()::Vector{Symbol} = [

#						:preserved_symmetries, 

#						:perturb_strength,  # only for plotting 

						:nr_perturb_strength,
						:max_perturb_strength,
						:nr_perturb_instances,
						]

function usedkeys(P::UODict)::Vector{Symbol} 

	@assert haskey(P, :preserved_symmetries)

	@assert haskey(P, :perturb_strength)


	uk = setdiff!(usedkeys(),[:perturb_strength])
	
	all_symms_preserved(P) && setdiff!(uk, [
								:nr_perturb_strength, # don't disregard!!
								:max_perturb_strength,
								:nr_perturb_instances
								]
					 )

	return uk

end 

obs_batch_1 = ["D110", "D111", "D113", "WannierGap"]
obs_batch_2 = ["D30", "D48","D123", "D125", "D127.1", "D127"]  

calc_observables = vcat(obs_batch_1, obs_batch_2)

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function perturb_strengths(P::UODict)::Vector{Float64}

	# the first must be zero 

	LinRange(0, P[:max_perturb_strength], P[:nr_perturb_strength]) 

#	all_symms_preserved(P) && return zeros(P[:nr_perturb_strength])


end  

function nr_perturb_instances(P::UODict)::Int 

	P[:nr_perturb_instances]

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

function get_perturb_on_mesh(P::UODict, seed::Int...;
														 kwargs...
														 )::Vector{Array{ComplexF64,4}} 
	
	get_perturb_on_mesh(preserved_symmetries(P), 
											nr_kPoints(P),
											kPoint_start(P),
											nr_perturb_instances(P),
											seed...;
											kwargs...)

end 

function get_perturb_on_mesh(
										 symms::AbstractString, 
										 args...; kwargs...
										 )::Vector{Array{ComplexF64,4}} 

	all_symms_preserved(symms) && return [zeros(ComplexF64,4,4,n-1,n-1)] 

	return get_perturb_on_mesh_(symms, args...; kwargs...)

end  

function get_perturb_on_mesh_(
										 symms::AbstractString, 
										 n::Int, k0::Real,
										 nr_pert::Int,
										 seed::Int
										 )::Vector{Array{ComplexF64,4}} 

	Random.seed!(seed)  

	return get_perturb_on_mesh_(symms, n, k0, nr_pert)

end  

function get_perturb_on_mesh_(
										 symms::AbstractString, 
										 n::Int, k0::Real,
										 nr_pert::Int,
										 )::Vector{Array{ComplexF64,4}}

	@assert !all_symms_preserved(symms) # safety -- methods changed
	
	randperts = [CalcWLO.get_perturb_on_mesh_(symms,n,k0) for i=1:nr_pert] 

	for i=2:nr_pert 

		@assert !isapprox(randperts[1][1], randperts[i][1], atol=1e-8) 

	end 

	return randperts 
 
end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function get_data_args(psiH::AbstractArray{<:Number,4},
											 results::AbstractDict,
											 f::Function=identity
											 )::Vector{Tuple{Array,Bool,Int,Bool}}

	r = any(in(keys(results)), obs_batch_2)

	return [(f(psiH), true, 1, r), (f(psiH), false, 1, r), 
					(f(psiH), true, 2, r), (f(psiH), false, 2, r),
					]

end  

function get_data_dir1(data::AbstractVector, dir1::Int
											)::AbstractVector

	view(data, 2dir1-1:2dir1)

end 


function get_data(psiH::AbstractArray{ComplexF64,4}, 
									results::AbstractDict;
									parallel::Bool=false
									)::Vector 
	
	if parallel

		args = get_data_args(psiH, results, copy)

		return pmap(Base.splat(WLO.get_wlo_data_mesh), args)

	else 

		args = get_data_args(psiH, results)

		return map(Base.splat(WLO.get_wlo_data_mesh), args)

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
	
	s = "_legend" 

	return length(obs)>length(s) && obs[end-length(s)+1:end]==s 

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


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function get_label(c::CartesianIndex, args...)::String

	get_label(c.I, args...)

end 

function get_label(I::Tuple{Vararg{Int}}, legend::OrderedDict)::String

	get_label(I, collect(values(legend)))

end 


function get_label(I::Tuple{Vararg{Int}}, 
									 legend::AbstractVector{<:AbstractVector})::String

	#	[string(leg,": \$",ylab[i],"\$") for (i,(leg,ylab)) in zip(c.I,pairs(legend))]

	#	return "\$"*join(map(Base.splat(getindex), zip(legend, I)),",")*"\$" 
	
	isempty(I) ? "" : "\$"*join(getindex.(legend, I),",")*"\$" 

end 

function obs_CI(Y::AbstractArray{<:Real,N}, 
								obs_i::Int, ax::AbstractVector{Int}=2:N 
							 )::CartesianIndex where N
		
	CartesianIndices(axes(Y)[ax])[min(end,obs_i)]

end 


function extract_data_and_lab(data::AbstractDict, obs::AbstractString,
													args...; kwargs...)::Tuple

	extract_data_and_lab(data[fn_nolegend(obs)], data[fn_legend(obs)], 
											 args...; kwargs...)

end 


function extract_data_and_lab(Y::AbstractArray{<:Number},
															legend::OrderedDict,
															obs_i::Int;
															kwargs...
															)::Tuple{Vector{<:Number}, String } 

	extract_data_and_lab(Y, collect(values(legend)), obs_i; kwargs...)

end 


function extract_data_and_lab(Y::AbstractArray{<:Number},
															axvars::AbstractVector{<:AbstractVector{<:AbstractString}},
															obs_i::Int; 
															transf_data::Function=copy,
															kwargs...
															)::Tuple{Vector{<:Number}, String } 
	

	c = obs_CI(Y, obs_i)

	return transf_data(view(Y, :, c)), get_label(c, axvars)

end  

function extract_data_and_lab(Y::AbstractArray{<:Number},
															legend::OrderedDict,
															obs_gr::Union{AbstractString,
																								<:AbstractVector{<:AbstractString}
																								},
															args...;
															kwargs...
															)::Tuple{Vector{Vector{<:Number}},
																					Vector{String},
																					String,
																					} 

	extract_data_and_lab(Y, 
											 (collect(keys(legend)), collect(values(legend))),
											 obs_gr,
											 args...;
											 kwargs...)

end 



function extract_data_and_lab(Y::AbstractArray{<:Number},
															legend::Tuple,
															obs_group::AbstractString,
															args...; kwargs...
															)::Tuple{Vector{Vector{<:Number}},
																					Vector{String},
																					String,
																					} 


	extract_data_and_lab(Y, legend, split_groups(obs_group), args...; kwargs...)

end 

function closest_periodic_iB(a::Real,
														 B::AbstractArray{<:Real},
														 T::Real)::Int

	D = Float64[0.0, 0.0]  # [min_dist, current_dist] 

	imin = 1 

	Utils.dist_periodic_!(D, 1, a, B[imin], T)

	for i=2:lastindex(B)

		Utils.dist_periodic_!(D, 2, a, B[i], T)

		D[2]<D[1] || continue 
	
		imin = i

		D[1] = D[2] 

	end 

	return imin

end 



function closest_periodic_iB(A::AbstractArray{<:Real},
														B::AbstractArray{<:Real},
														T::Real
														)::Int
	
	bcount = zeros(Int,length(B))

	for a in A 

		bcount[closest_periodic_iB(a,B,T)] += 1 

	end 

	return argmax(bcount)

end 
function closest_periodic_b(A::AbstractArray{<:Real},
														B::AbstractArray{<:Real},
														T::Real
														)::Float64

	B[closest_periodic_iB(A,B,T)]

end  




function extract_data_and_lab(Y::AbstractArray{<:Number,N},
														 (axnames,axvars)::Tuple{
												 <:AbstractVector{<:AbstractString},
												 <:AbstractVector{<:AbstractVector{<:AbstractString}}
												 },
															obs_groups::AbstractVector{<:AbstractString},
															obs_i::Int;
															transf_data::Function=copy!,
															#inplace!
															out_type::DataType=Float64,
															kwargs...
															)::Tuple{Vector{Vector{<:Number}},
																					Vector{String},
																					String,
																					} where N 


	I = Vector{Union{Base.OneTo{Int},Int}}(undef, N-1) 

	fullax = map(in(obs_groups), axnames)
	partax = findall(!,fullax)  

	c = obs_CI(Y, obs_i, partax .+ 1)
	
	I[fullax] .= axes(Y)[(2:N)[fullax]]
	I[partax] .= c.I

	fixedlabel = get_label(c, axvars[partax]) 

	nr_curves = mapreduce(length,*,I)

	chlabels = Vector{String}(undef,nr_curves)



	ys = [Vector{out_type}(undef,size(Y,1)) for i=1:nr_curves] 


	for (n,i) in enumerate(Iterators.product(I...)) 

		y0 = closest_periodic_b(view(Y,:,i...),[0,0.5],1)

		chlabels[n] = myPlots.join_label(get_label(i[fullax], axvars[fullax]),
																		 string(y0))

#		chlabels[n] = get_label(i[fullax], axvars[fullax])  


		if y0≈0 
			transf_data(ys[n], view(Y,:,i...))
		else 
			transf_data(ys[n], Utils.closest_periodic_shifted_a.(view(Y,:,i...).-y0,0,1))

		end 

		y0count = 0 

	end 

	return ys, chlabels, fixedlabel

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




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



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

join_groups =  Base.Fix2(join," ")
split_groups = Base.Fix2(split," ")

function combs_groups(g::AbstractVector{<:AbstractString
																			 }=["dir1","sector","stat"]
										 )::Vector{String}

	map(join_groups, Combinatorics.powerset(g,1))

end 


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

	for eq in ("D123","D127.1","D48","D113","D127")

		in(eq,obs) || continue 

		results[eq] = zeros(L,2,2,2)

		add_legend!(results, eq, ["dir1","sector","stat"])

	end 


	for eq in ("D110","D111","D125","D30",)#"D48")

		in(eq,obs) || continue 

		results[eq] = zeros(L,2,2)

		add_legend!(results, eq, ["dir1", "stat"])

	end 

#	for eq in ("D30","D48",) 
#
#		in(eq,obs) || continue 
#		
#		results[eq] = zeros(L,2)
#
#		add_legend!(results, eq, ["stat"])
#
#	end 

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

		obs in ("D123","D48")  && return ["+", "-"]

		if obs in ("D113","D127", "D127.1") 
			
			return ["\\mathrm{occup}","\\mathrm{unocc}"]

			# tex command \text{...} not understood in matplotlib 
			
		end 

	end  

	error("Wrong input")

end 

function set_result_one!(results::AbstractDict{String,Any},
												 obs::AbstractString,
												 numerical_result::T,
												 trial::Int, i_ps::Int, inds...
												 )::T where T<:AbstractVector{<:Real} 


	
	I = map(inds,values(results[fn_legend(obs)])) do val,vals 

				val isa Int && return val 

				i = findfirst(==(string(val)),vals)

				if isnothing(i) 

					i = findfirst((occursin(val,V) for V=vals))

				end 

				@assert i isa Int  

				return i 

			end 
				
	SignalProcessing.run_m1!(results[obs], trial, numerical_result, 
													 i_ps, I..., :)

	return numerical_result 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function set_results!(results::AbstractDict,
													nk::Int, k0::Real,
													trial::Int, i_ps::Int,
													data)::Nothing

	for dir1=1:2 
					
		set_results_onedir!(results, nk, k0, trial, i_ps, dir1, 
												get_data_dir1(data, dir1)...)

	end  

	set_results_twodir!(results, nk, k0, trial, i_ps, 
											get_data_dir1(data, 1)..., 
											get_data_dir1(data, 2)... 
											)

end   


function set_results_twodir!(results::AbstractDict, nk::Int, k0::Real,
													trial::Int,i_ps::Int,
													(eigW1_occup_x,nus2pm_xy)::Tuple{
																<:AbstractVector{<:AbstractArray},
																<:AbstractVector{<:AbstractArray{<:Real,3}}
																			},
													(eigW1_unocc_x,eta2pm_xy)::Tuple{
																<:AbstractVector{<:AbstractArray},
																<:AbstractVector{<:AbstractArray{<:Real,3}}
																			},
													(eigW1_occup_y,nus2pm_yx)::Tuple{
																<:AbstractVector{<:AbstractArray},
																<:AbstractVector{<:AbstractArray{<:Real,3}}
																			},
													(eigW1_unocc_y,eta2pm_yx)::Tuple{
																<:AbstractVector{<:AbstractArray},
																<:AbstractVector{<:AbstractArray{<:Real,3}}
																			},

													)::Nothing  


	if haskey(results,"D127.1") 

	#	D127.1: polariz = sum over Wannier sectors  
	
		p1 = WLO.polariz_fromSubspaces1(eigW1_occup_x)
		p1_ = WLO.add_nupm(nus2pm_yx)

		set_result_one!(results, "D127.1", 
										abs.(WLO.wcc_stat_diff!(p1_,p1,[0,0.5])), 
										trial, i_ps, "x", "occup")


		WLO.polariz_fromSubspaces1!(p1, eigW1_unocc_x)  
		WLO.add_nupm!(p1_, eta2pm_yx)
	
		set_result_one!(results, "D127.1", 
										abs.(WLO.wcc_stat_diff!(p1_,p1,[0,0.5])), 
										trial, i_ps, "x", "unocc")


		WLO.polariz_fromSubspaces1!(p1, eigW1_occup_y)
		WLO.add_nupm!(p1_, nus2pm_xy) 

		set_result_one!(results, "D127.1", 
										abs.(WLO.wcc_stat_diff!(p1_,p1,[0,0.5])), 
										trial, i_ps, "y", "occup")


		WLO.polariz_fromSubspaces1!(p1, eigW1_unocc_y)  
		WLO.add_nupm!(p1_, eta2pm_xy)
	
		set_result_one!(results, "D127.1", 
										abs.(WLO.wcc_stat_diff!(p1_,p1,[0,0.5])), 
										trial, i_ps, "y", "unocc")
	end 






	return  



end 

function set_results_onedir!(results::AbstractDict, nk::Int, k0::Real,
													trial::Int,i_ps::Int,dir1::Int,
													(eigW1_occup,nus2pm)::Tuple{
																<:AbstractVector{<:AbstractArray},
																<:AbstractVector{<:AbstractArray{<:Real,3}}
																			},
													(eigW1_unocc,eta2pm)::Tuple{
																<:AbstractVector{<:AbstractArray},
																<:AbstractVector{<:AbstractArray{<:Real,3}}
																			},
													)::Nothing

#	WannierGap 					
	if haskey(results,"WannierGap")

		gap = min(WLO.WannierGap_fromSubspaces(eigW1_occup),
							WLO.WannierGap_fromSubspaces(eigW1_unocc))
	
	
		SignalProcessing.run_m1!(results["WannierGap"], trial, gap, i_ps, dir1)

	end 



	# -------------------------------------------- #

	if any(in(keys(results)), obs_batch_1)

		p1occup = MODEL.polariz_fromSubspaces1(eigW1_occup,dir1) # vector
		p1unocc = MODEL.polariz_fromSubspaces1(eigW1_unocc,dir1) # vector 


		if haskey(results,"D110")
	
			#	D110: p1occup-p1unocc ==0
				SignalProcessing.run_m1!(results["D110"],trial,
										abs.(WLO.wcc_stat_diff(p1occup,p1unocc)),
										i_ps,dir1,:)
			
		end 
		
		if haskey(results,"D111")

		#	D111: p1occup+p1unocc ==0
			SignalProcessing.run_m1!(results["D111"],trial,
									abs.(WLO.wcc_stat_sum(p1occup,p1unocc)),
									i_ps,dir1,:)
		
		end 	
	
		if haskey(results,"D113")

		#	D113: p1occup or p1unocc in [0,0.5] 

			set_result_one!(results, "D113",
											abs.(WLO.wcc_stat(p1occup,[0,0.5])),
											trial, i_ps, dir1, "occup") 

			set_result_one!(results, "D113",
											abs.(WLO.wcc_stat(p1unocc,[0,0.5])),
											trial, i_ps, dir1, "unocc") 

		end 	
	
		if haskey(results,"D125")||haskey(results,"D127") 


			p1occup_ = MODEL.add_nupm(nus2pm,3-dir1)
			p1unocc_ = MODEL.add_nupm(eta2pm,3-dir1)


			if haskey(results,"D127")
				# occupied polarization quantized. Should be identical to D113
		
				set_result_one!(results, "D127",
												abs.(WLO.wcc_stat(p1occup_, [0,0.5])),
												trial, i_ps, dir1, "occup") 


				set_result_one!(results, "D127",
												abs.(WLO.wcc_stat(p1unocc_, [0,0.5])),
												trial, i_ps, dir1, "unocc")

				
			end 

			if haskey(results,"D125")
				#	D125: total polarization zero  
				# Should be identical to D111
			
					SignalProcessing.run_m1!(results["D125"], trial,
											abs.(WLO.wcc_stat_sum(p1occup_,p1unocc_,[0,0.5])),
										 i_ps,dir1,:)
				
			end 
		end 

	end 


	# -------------------------------------------- #
	if haskey(results,"D123")  

		for sector in 1:2 # plus and minus  sectors 

#	D123: nu2pm == eta2pm 
			
			set_result_one!(results, "D123",
											abs.(WLO.wcc_stat_diff(nus2pm[sector], eta2pm[sector])), 
											trial, i_ps, dir1, sector)

		end # sector 

	end



#		D30 for (d1,d2)==(1,2), one Wannier band per sector

	if haskey(results,"D30")

			nu2_minus_mk = selectdim(nus2pm[2],dir1+1,
								WLO.ind_minusk.(axes(nus2pm[2],dir1+1),nk,k0)
								)

			SignalProcessing.run_m1!(results["D30"], trial,
									abs.(WLO.wcc_stat_diff(nus2pm[1],nu2_minus_mk)),
								 i_ps,dir1,:)

	end 

#	D48 for dir1==2 needed. Calc also dir1==1
		
	if haskey(results,"D48")

		for sector in 1:2 

#nu_y recoreded at dir1==1=x
#nu_x recoreded at dir1==2=y

			SignalProcessing.run_m1!(results["D48"], trial,
									abs.(WLO.wcc_stat(view(nus2pm[sector],:),
																		 [0,0.5])),
									i_ps, dir1, sector, : )


		end 

#
	end 
	

	return 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


prep_obs(::Nothing=nothing)::Vector{String} = String[]

function prep_obs(obs::AbstractString)::Vector{String} 
	
	fn_nolegend(obs) in calc_observables || return prep_obs() 

	return vcat(fn_nolegend(obs),fn_legend(obs),xxlabel())

end 

function prep_obs(obs::AbstractVector{<:AbstractString})::Vector{String} 

	mapreduce(prep_obs, vcat, obs; init=prep_obs())

end 

function prep_obs(args...)::Vector{String} 
	
	mapreduce(prep_obs, vcat, args; init=prep_obs())

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

get_target = prep_obs ∘ Helpers.f_get_target(:observables)	   

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
	parallel=false#nprocs()>=4
	parallel2=nprocs()>=4

	strengths = perturb_strengths(P)

	zero_strength, rest_strengths = Base.Iterators.peel(strengths)
	s1, s2e = Base.Iterators.peel(eachindex(strengths))

	@assert iszero(zero_strength)


	results = init_results(strengths, get_target(target; kwargs...))


#	return Dict{String,Any}(k=>v for (k,v)=pairs(results) if isauxfile(k)) 


	all(isauxfile, keys(results)) && return results



	# ------ no perturbation --------- #

	set_results!(results, nk, k0, 1, s1,  
										 get_data(MODEL.get_psiH(P, nk, k0, 
#																						 perturb2, #### test 
																						 ), 
															results;
															parallel=parallel2))


	if all_symms_preserved(P) 

		for (k,v) in pairs(results)

			isauxfile(k) && continue
	
			for s in s2e 

				copy!(selectdim(results[k],1,s), selectdim(v,1,s1))

			end 

		end 

		return results

	end 	

	# ------------- with perturbation -------------- # 

	trials = get_perturb_on_mesh(P)#, 3268) # seed ensures results recovered

	
	if parallel #---- parallel evaluation ---#

		#error() #################### test 

		data_all = pmap(Iterators.product(trials,rest_strengths)) do pert

			get_data(MODEL.get_psiH(P, nk, k0, pert...), results)

		end # get_data contains 4 jobs  


		for (I,d) in zip(Iterators.product(eachindex(trials),s2e), data_all)

			set_results!(results, nk, k0, I..., d)

		end  


	else  #--------#

		for (j,perturb0) in enumerate(trials) # don't parallelize this j-loop! 
	
			for (i,ps) in zip(s2e,rest_strengths) 
	
				data = get_data(MODEL.get_psiH(P, nk, k0, 
																			 perturb0, ps, 
#																			perturb0*ps + perturb2 ####### test 
																			), results;
												parallel=parallel2)
	
				set_results!(results, nk, k0, j, i, data)
	
			end  
	
		end # ---------#  

	end 

	return results

end 








#===========================================================================#
#
#
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
end  # module ChecksWLO

