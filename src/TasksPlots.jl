module TasksPlots
#############################################################################


#import LinearAlgebra,QuadGK,Optim

#import myLibs: Utils, ComputeTasks, Algebra, Parameters, Lattices
#import myLibs.Parameters:UODict 




import myLibs: Utils, ComputeTasks, SignalProcessing 

using myLibs.ComputeTasks: CompTask   

import myPlots
using myPlots: PlotTask 

import ..ChecksWLO, ..Helpers, ..CalcWLO, ..WLO, ..RibbonWLO

using ..Helpers.hParameters: Calculation  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function apply_log10abs(src::AbstractVector{<:Real})::Vector{Float64}

	apply_log10abs!(Vector{Float64}(undef, length(src)), src)

end 

function apply_log10abs!(dst::AbstractVector{Float64},
											src::AbstractVector{<:Real}
											)::AbstractVector{Float64}

	m::Float64 = -18

	nz = !Utils.fSame(0, 10^m)


	i = Union{Nothing,Int}[0,0]

	i[1] = findfirst(nz, src)

	isnothing(i[1]) && return fill!(dst, m) 


	fill!(view(dst,firstindex(dst):i[1]), log10(abs(src[i[1]])))

	while !isnothing(setindex!(i, findnext(nz, src, i[1]+1), 2)[2])

		map!(SignalProcessing.LinearInterp1D(dst[i[1]], log10(abs(src[i[2]])), 
																				 i[1], i[2]),
				 view(dst, i[1]+1:i[2]), i[1]+1:i[2])

		i[1] = i[2] 

	end 


	fill!(view(dst, i[1]+1:lastindex(dst)), dst[i[1]])


	m1 = minimum(dst)
	
	m1<m && @warn "Update the minimum: value $m1 encountered"

	return dst 

end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function init_sliders_obs(observables::AbstractVector{<:AbstractString}
												 )::Vector
													
	[ ("observables", observables), ("obs_index", 8), ]


end 
#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_ylim(::Nothing, ::Nothing, ydata::AbstractArray, args...
									)::Vector{Float64}

									
	Utils.extend_limits(minimum(minimum,ydata), maximum(maximum,ydata), args...)

end 

function get_ylim(::Nothing, given_ymax::Real, ydata::AbstractArray, args...
									)::Vector{Float64}

	ylim = Utils.extend_limits(minimum(minimum,ydata), given_ymax, args...)

	ylim[2] = given_ymax 

	return ylim 

end 
function get_ylim(given_ymin::Real, ::Nothing, ydata::AbstractArray, args...
									)::Vector{Float64}

	ylim = Utils.extend_limits(given_ymin, maximum(maximum, ydata), args...)

	ylim[1] = given_ymin 

	return ylim 

end 

function get_ylim(given_ymin::Real, given_ymax::Real, args...
									)::Vector{Float64}

	[given_ymin, given_ymax]

end 

function get_ylim(P::AbstractDict, args...)::Vector{Float64}

	get_ylim(get(P,"ymin",nothing), get(P,"ymax",nothing), args...) 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function f_extract_data_and_lab(P::AbstractDict; kwargs...)::Function

	function apply_rightaway(data::AbstractDict, args... 
													 )::Tuple{Tuple,Tuple}

		(ChecksWLO.xxlabel(data),
		 ChecksWLO.extract_data_and_lab(data, P["obs"], 
																		get(P,"obs_group","-"), 
																		get(P,"obs_i",1);
																		kwargs...)
		 )

	end 

end 

function f_extract_data_and_lab_atPS(P::AbstractDict; 
																		 atol::Float64=1e-10,
																		 kwargs...)::Function

	PS = get(P, "perturb_strength", 0)

	function F2_(((x,xlabel),(ys,chlabels,fixedlab)
								)::Tuple{Tuple,Tuple}
							 )::Tuple

		if isapprox(extrema(x)..., atol=atol)

			for y in ys 
				if !isapprox(extrema(y)...,atol=atol)
					@warn "x identical but y different"
					break
				end 
			end 

			return (PS, xlabel, first.(ys), chlabels, fixedlab)

		else 

			y_ps = [SignalProcessing.Interp1D(x,y,1,PS) for y=ys] 

			if any(isnan, y_ps) || any(isinf, y_ps) 

				@show ys y_ps x PS 

				error("NaN/Inf found in interpolated ys")

			end 

			return (PS, xlabel, y_ps, chlabels, fixedlab)

		end 

	end 

	return F2_ ∘ f_extract_data_and_lab(P; kwargs...)

end 









#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




tex(lab::AbstractString)::String = "\$$lab\$"   

function plotdict_WCC(P::AbstractDict, 
											task::CompTask,
											obs::AbstractString,
											d::Int; 
											)::Dict{String,Any}

	data = task.get_data(P; fromPlot=true, mute=false, target=obs)  

	ks,klabels = CalcWLO.xxlabel(data)

	legend = data[CalcWLO.fn_legend(obs)] 

	y0 = transpose(selectdim(data[obs],2,d))

	#yave = WLO.wcc_stat!(sum(eachrow(y0)),[0,0.5])[1] 

	yave = y0[2]

#	@show yave 
	ks ./= pi 

	return Dict{String,Any}(

					"xs"=> repeat(transpose(ks),outer=(size(y0,1),3)),

					"ys"=> hcat(y0 .- 1, y0, y0 .+ 1),
					
					"labels" => map(tex, legend["sector"]),

					"xlim"=> Utils.extend_limits(extrema(ks),1e-3),

					"ylim"=> (-0.7, 0.7),

					"xlabel" => tex(klabels[d]*"/\\pi"),

					"ylabel" => tex(legend["dir"][d]),
		
					"yline" => Utils.closest_periodic_b(yave,[0,0.5],1),

					)

end 

function plotdict_checkzero(
												 (x,xlabel)::Tuple{<:AbstractVector,<:AbstractString},
												 (ys,chlabels,fixedlab
													)::Tuple{<:AbstractVector{<:AbstractVector},
																	 <:AbstractVector{<:String},
																	 <:AbstractString},
												 P::AbstractDict;
												 transf_lab::AbstractString="",
												 kwargs...
												 )::Dict{String,Any} 


	out = Dict{String,Any}(

		"xs"=> [x for y=ys],

		"ys" => ys, 

		"xlim" => extrema(x),

		"ylim" => get_ylim(P, ys),

		"labels" => chlabels,

		"ylabel" => myPlots.join_label(transf_lab*P["obs"], fixedlab),

		"xlabel" => xlabel,

		)

	ComputeTasks.add_line!(out, P, "perturb_strength", "x")

	return out 

end  


function plot_check_zero(P::AbstractDict, 
												 task::CompTask; 
												 kwargs...)::Dict{String,Any} 

	F = f_extract_data_and_lab(P; kwargs...)

	data = F(task.get_data(P; fromPlot=true, mute=false, target=P["obs"]))

	return plotdict_checkzero(data..., P; kwargs...)

end 

function plot_check_zero_multi(P::AbstractDict, 
															 task::CompTask, 
															 out_dict::AbstractDict; 
															 kwargs...)::Dict{String,Any} 

	F = f_extract_data_and_lab(P; kwargs...)


	data = task.get_data(P; fromPlot=true, target=P["obs"], apply_rightaway=F)
	
	curves = vcat((d[2][1] for d in data)...)

	chlabs = vcat((myPlots.join_label.(d[2][2], [string(y)]) for (d,y)=zip(data,out_dict["y"]))...)

	fixedlab = myPlots.join_label(data[1][2][3], out_dict["ylabel"])


	return plotdict_checkzero(data[1][1],
														(curves,chlabs,fixedlab),
														P; 
														kwargs...)
end  


function add_default_obs(observables::AbstractVector{<:AbstractString}
												 )::Function 

	P0 = Dict{String,Any}("obs"=>first(observables))

	return function adjust_params(P::AbstractDict)::AbstractDict

		(haskey(P,"obs") && in(P["obs"],observables)) ? P : merge(P,P0)

	end 

end  


function get_dir1(P::AbstractDict)::Int 

	min(2, get(P,"obs_i",1))  

end 

function get_dir2(P::AbstractDict)::Int 

	3-get_dir1(P)

end 


function nr2str(x::Real,n::Int=3; fixedsize::Bool=true)::String 

	abs(x)<1e-10 && return "0.0"

	d = Int(n-ceil(log10(abs(x))))

	s = string(round(x, digits=d))[1:min(end,d+n)]

	return fixedsize ? rpad(s,d+n-1,'0') : s

end 
	

function Wannier_gap_text(gap::Real)::String  

#	gap<1e-10 && return "(gap\$=0.0\$)"

	string("(gap\$=", nr2str(gap, 3)[1:min(end,5)], "\$)")

end 

function wcc_stat_text(stat::AbstractVector{<:Real})::String  

	wcc_stat_text(nr2str.(stat, 3))

end  

function wcc_stat_text(stat::AbstractVector{<:AbstractString})::String  

	tex("("* join(stat,",")*")")

end 


#===========================================================================#
#
function CheckZero(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#


	observables_ = intersect(observables, ChecksWLO.calc_observables)



	task = CompTask(Calculation("Check zeros",
															ChecksWLO, init_dict; 
															observables=observables_, kwargs...))


	function plot(P::AbstractDict)::Dict{String,Any} 

		plot_check_zero(P, task; 
										transf_lab="\$\\log_{10}\$",
										transf_data=apply_log10abs!
										)

	end

	return PlotTask(task, 
									[init_sliders_obs(observables_);
									 [("obs_group", ChecksWLO.combs_groups()),	("ylim",)]], 
									"Curves_yofx", plot∘add_default_obs(observables_))
end 

#===========================================================================#
#
function RibbonWannierBands1(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	obs = "WannierBands1"

	observables_ = filter(==(obs), observables)

	task = CompTask(Calculation("Wannier bands ribbon",
															RibbonWLO, init_dict; 
															observables=observables_, kwargs...))


	function plot(P::AbstractDict)::Dict{String,Any} 

		d = get_dir1(P) 

		data = task.get_data(P; fromPlot=true, mute=false, target=obs)  
	
		x,xlabel = RibbonWLO.xxlabel(data)
	
		legend = data[RibbonWLO.fn_legend(obs)] 
	
		y0 = selectdim(data[obs],2,d) #Vector
	
		gap = 1.0 

		for i = 1:length(y0)-1 
			
			gap = min(gap,
								Utils.reduce_dist_periodic(min, y0[i], view(y0, i+1:length(y0)), 1)
								)

		end 

		xlim = Utils.extend_limits(extrema(x),1e-3)

		xlim[2] = max(xlim[2],xlim[1]+1)

		out = Dict{String,Any}(
	
					 "xs" => [x,x,x],
	
					 "ys" => [y0 .- 1, copy(y0), y0 .+ 1],
						
					 "labels" => ["-1","0","+1"],
	
						"xlim"=> xlim,
	
						"ylim"=> (-0.7, 0.7),
	
						"xlabel" => xlabel,
	
						"ylabel" => tex(legend["dir"][d]) * " " * Wannier_gap_text(gap),

						"yline" => 0.5,
	
						)
	
	
	
			return out 
	
		end 

	return PlotTask(task, ("obs_index", 2), "Scatter", plot,)
end 

#===========================================================================#
#
function WannierBands1(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	obs = "WannierBands1"

	observables_ = filter(==(obs), observables)

	task = CompTask(Calculation("Wannier bands",
															CalcWLO, init_dict; 
															observables=observables_, kwargs...))


	function plot(P::AbstractDict)::Dict{String,Any} 

		out = plotdict_WCC(P,task,obs,get_dir1(P))

		gap = Utils.reduce_dist_periodic(min, eachrow(out["ys"])..., 1)

		out["ylabel"] *= " "*Wannier_gap_text(gap)

		return out 

	end 


	return PlotTask(task, ("obs_index", 2), "Scatter", plot,)
end 

#===========================================================================#
#
function WannierBands2(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	obs = "WannierBands2"

	observables_ = filter(==(obs), observables)

	task = CompTask(Calculation("Nested Wannier bands",
															CalcWLO, init_dict; 
															observables=observables_, kwargs...))


	function plot(P::AbstractDict)::Dict{String,Any} 

		out = plotdict_WCC(P,task,obs,get_dir2(P))

		for i=1:2 

			stat = -apply_log10abs(WLO.wcc_stat(selectdim(out["ys"],1,i), [0,0.5]))

			out["labels"][i] *= " " * wcc_stat_text(stat)

		end  




		out["ylabel"] *= ",  L:\$-\\log_{10}\$" * wcc_stat_text(["\\mu","\\sigma"])
	
		out["ylabel"] *= "="*wcc_stat_text(-apply_log10abs(WLO.wcc_stat!(sum(eachrow(out["ys"])),[0,0.5])))

		return out 

	end 


	return PlotTask(task, ("obs_index",2), "Scatter", plot,)
end 

#===========================================================================#
#
function WannierGap(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	obs = "WannierGap"

	observables_ = filter(==(obs), observables)

	task = CompTask(Calculation(obs, ChecksWLO, init_dict; 
															observables=observables_,
															kwargs...))

	function plot(P::AbstractDict)::Dict{String,Any} 

		out = plot_check_zero(merge(P, Dict("obs"=>obs)), task)
		
		return setindex!(out, [0,0.07], "ylim")

	end

	return PlotTask(task, 
									[init_sliders_obs(observables_);
									 ("obs_group", ChecksWLO.combs_groups())
									 ], "Curves_yofx", plot)

end 



#===========================================================================#
#
function WannierGap_atYs(init_dict::AbstractDict;
											 Y::Symbol,
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	obs = "WannierGap"
	
	observables_ = filter(==(obs), observables)

	C = Calculation("$obs for several $Y", ChecksWLO, init_dict; 
									observables=observables_, kwargs...)

	task, out_dict, = ComputeTasks.init_multitask(C, [Y=>1], [1=>""])

	P0 = Dict("obs"=>obs, "obs_group"=>"-")

	function plot(P::AbstractDict)::Dict{String,Any} 
	
		out = plot_check_zero_multi(merge(P, P0), task, out_dict)

		return setindex!(out, [0,0.07], "ylim")
		
	end

	return PlotTask(task, init_sliders_obs(observables_), "Curves_yofx", plot)

end 


#===========================================================================#
#
function CheckZero_atYs(init_dict::AbstractDict;
											 Y::Symbol,
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	observables_ = intersect(observables, ChecksWLO.calc_observables) 

	C = Calculation("Check zeros for several $Y", ChecksWLO, init_dict; 
									observables=observables_, kwargs...)

	task, out_dict, = ComputeTasks.init_multitask(C, [Y=>1], [1=>""])

	P0 = Dict{String,Any}("obs_group"=>"-") 

#		plot∘add_default_obs(observables_)∘Base.Fix1(merge,P0) 

	function plot(P::AbstractDict)::Dict{String,Any} 

		P1 = if haskey(P,"obs") && in(P["obs"],observables_)

				copy(P)

			else 

				merge(P,Dict{String,Any}("obs"=>first(observables)))

			end  


		return plot_check_zero_multi(merge!(P1,P0), task, out_dict;
										transf_lab="\$\\log_{10}\$",
										transf_data=apply_log10abs!
										)

	end

	return PlotTask(task, 
									[init_sliders_obs(observables_); ("ylim",)], 
									"Curves_yofx", 
									plot)

end 







#===========================================================================#
#
function CheckZero_atPS_vsX(init_dict::AbstractDict;
											 X::Symbol,
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	observables_ = intersect(observables, ChecksWLO.calc_observables) 

	C = Calculation("Check zeros at fixed PS vs. $X", 
									ChecksWLO, init_dict; 
									observables=observables_, kwargs...)

	task, out_dict, = ComputeTasks.init_multitask(C, [X=>1])#, [1=>""])

	if X==:nr_kPoints 
		out_dict["x"] = 1.0 ./ out_dict["x"]
		out_dict["xlabel"] = "1/("*out_dict["xlabel"]*")"
		out_dict["xlim"] = [0,maximum(out_dict["x"])]
	else 
		out_dict["xlim"] = extrema(out_dict["x"])
	end 

	P1 = Dict{String,Any}("obs"=>first(observables_))

	#d3(d::Tuple)::Float64 = only(d[3]) 

	function plot(P2::AbstractDict)::Dict{String,Any} 

		P = in(get(P2,"obs",0), observables_) ? P2 : merge(P2,P1)

		transf_lab="\$\\log_{10}\$" 

		kwargs = (transf_lab="\$\\log_{10}\$",
							transf_data=apply_log10abs!
							)

		F = f_extract_data_and_lab_atPS(P; kwargs...)

		data = task.get_data(P; fromPlot=true, target=P["obs"], apply_rightaway=F)
		
		chlabs = data[1][4] 

		ys = [zeros(length(data)) for l=chlabs]

		for (i,d) in enumerate(data)
			for (j,y) in enumerate(d[3])

				ys[j][i]=y

			end 
		end 


#		ys = mapreduce(Base.Fix2(getindex,3),hcat,data)


		ylabel = myPlots.join_label(transf_lab*P["obs"],
																data[1][5], 
																"PS="*nr2str(data[1][1],2),
																						 )

		out = Dict{String,Any}(

							"xs" => [out_dict["x"] for y=ys],
							
							"ys" => ys,

							"xlim" => out_dict["xlim"],

							"ylim" => get_ylim(P, ys),


							"labels" => chlabs,

							"ylabel" => ylabel,

							"xlabel" => out_dict["xlabel"],


						)


#										transf_lab="\$\\log_{10}\$",
#										transf_data=apply_log10abs!
#										)

		ComputeTasks.add_line!(out, P, X, "x") # also in plotdict_checkzero

		return out 
	end

	return PlotTask(task,  
									[init_sliders_obs(observables_);
									 [("obs_group", ChecksWLO.combs_groups()),	("ylim",)]], 
									"Curves_yofx", 
									plot)

end 


#===========================================================================#
#
function CheckZero_atPS_vsX_atYs(init_dict::AbstractDict;
											 X::Symbol,Y::Symbol,
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	observables_ = intersect(observables, ChecksWLO.calc_observables) 

	C = Calculation("Check zeros at fixed PS for several $Y",
									ChecksWLO, init_dict; 
									observables=observables_, kwargs...)

	task, out_dict, construct_Z, = ComputeTasks.init_multitask(C, [X=>1, Y=>1])

	if X==:nr_kPoints 
		out_dict["x"] = 1.0 ./ out_dict["x"]
		out_dict["xlabel"] = "1/("*out_dict["xlabel"]*")"
		out_dict["xlim"] = [0,maximum(out_dict["x"])]
	else 
		out_dict["xlim"] = extrema(out_dict["x"])
	end 

	P0 = Dict{String,Any}("obs_group"=>"-") 

	P1 = Dict{String,Any}("obs"=>first(observables_))
		
	d3(d::Tuple)::Float64 = only(d[3])

	function plot(P2::AbstractDict)::Dict{String,Any} 

		P = in(get(P2,"obs",0),observables_) ? merge(P2,P0) : merge(P2,P1,P0)

		transf_lab="\$\\log_{10}\$" 

		kwargs = (transf_lab="\$\\log_{10}\$",
							transf_data=apply_log10abs!
							)

		F = f_extract_data_and_lab_atPS(P; kwargs...)

		data = task.get_data(P; fromPlot=true, target=P["obs"], apply_rightaway=F)
	
		for d in data
			@assert length(d[3])==length(d[4])==1
#			@assert isempty(only(d[4])) 
			@assert !isnothing(tryparse(Float64, only(d[4])))
		end 

		chlabels = [only(d[4]) for d in data] 

		labels = string.(out_dict["y"]) 
		
		obs = P["obs"]

		if length(chlabels)>1 && length(unique(chlabels))>1 
			
			for (i,l) in enumerate(labels)

				labels[i] = myPlots.join_label(l,chlabels[i])

			end 

		elseif !(tryparse(Float64,chlabels[1])≈0)

			obs = string("($obs-",chlabels[1],")")

		end 




		ys = collect.(eachcol(construct_Z(d3, data)["z"]))


		ylabel = myPlots.join_label(transf_lab*obs,#P["obs"],
																data[1][5], 
																"PS="*nr2str(data[1][1],2),
																						 )


		out = Dict{String,Any}(

							"xs" => [out_dict["x"] for y=ys],

							"xlabel" => out_dict["xlabel"],

							"ys" => ys,

							"xlim" => out_dict["xlim"],

							"ylim" => get_ylim(P, ys),

							"labels" => out_dict["y"],
		
							"ylabel" => ylabel,

							)

	ComputeTasks.add_line!(out, P, X, "x") # also in plotdict_checkzero

	return out 

#										transf_lab="\$\\log_{10}\$",
#										transf_data=apply_log10abs!
#										)

	end

	return PlotTask(task, 
									[init_sliders_obs(observables_); ("ylim",)], 
									"Curves_yofx", 
									plot)

end 


































































































































































































































































































#############################################################################
end # module TasksPlots 
