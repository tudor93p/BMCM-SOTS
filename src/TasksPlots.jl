module TasksPlots
#############################################################################


#import LinearAlgebra,QuadGK,Optim

#import myLibs: Utils, ComputeTasks, Algebra, Parameters, Lattices
#import myLibs.Parameters:UODict 




import myLibs: Utils, ComputeTasks, SignalProcessing 

using myLibs.ComputeTasks: CompTask   

import myPlots
using myPlots: PlotTask 

import ..ChecksWLO, ..Helpers, ..CalcWLO

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

	ks ./= pi 

	return Dict{String,Any}(

					"xs"=> repeat(transpose(ks),outer=(size(y0,1),3)),

					"ys"=> hcat(y0 .- 1, y0, y0 .+ 1),
					
					"labels" => map(tex, legend["sector"]),

					"xlim"=> extrema(ks),

					"ylim"=> (-0.7, 0.7),

					"xlabel" => tex(klabels[d]*"/\\pi"),

					"ylabel" => tex(legend["dir"][d]),

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


function Wannier_gap_text(gap::Real)::String  

	gap<1e-10 && return " (gap\$=0.0\$)"

	gap = string(round(gap, digits=Int(ceil(-log10(gap)))+2))[1:min(end,6)]

#	gap = rpad(gap, 6, "0")

	return " (gap\$=$gap\$)"

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

		out["ylabel"] *= 	Wannier_gap_text(gap)

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

		plotdict_WCC(P,task,obs,get_dir2(P))

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
function WannierGap_atY(init_dict::AbstractDict;
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
function CheckZero_atY(init_dict::AbstractDict;
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










































































































































































































































































































#############################################################################
end # module TasksPlots 
