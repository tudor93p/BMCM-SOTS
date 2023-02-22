module TasksPlots
#############################################################################


#import LinearAlgebra,QuadGK,Optim

#import myLibs: Utils, ComputeTasks, Algebra, Parameters, Lattices
#import myLibs.Parameters:UODict 




import myLibs: Utils, ComputeTasks, SignalProcessing 

using myLibs.ComputeTasks: CompTask   

import myPlots
using myPlots: PlotTask 

import ..ChecksWLO, ..Helpers 

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
													
	[ ("observables", observables), 
	 ("obs_index", 8), 
	 ]

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
		 ChecksWLO.extract_data_and_lab(data, P["obs"], P["obs_group"], P["obs_i"];
																		kwargs...)
		 )

	end 

end 


function plotdict_checkzero(
												 (x,xlabel)::Tuple{<:AbstractVector,<:AbstractString},
												 (ys,chlabels,fixedlab
													)::Tuple{<:AbstractVector{<:AbstractVector},
																	 <:AbstractVector{<:String},
																	 <:AbstractString},
												 P::AbstractDict=Dict();
												 transf_lab::AbstractString="",
												 kwargs...
												 )::Dict{String,Any} 

	Dict{String,Any}(

		"xs"=> [x for y=ys],

		"ys" => ys, 

		"xlim" => extrema(x),

		"ylim" => get_ylim(P, ys),

		"labels" => chlabels,

		"ylabel" => myPlots.join_label(transf_lab*P["obs"], fixedlab),

		"xlabel" => xlabel,

		)

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



#===========================================================================#
#
function CheckZero(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	task = CompTask(Calculation("Check zeros",
															ChecksWLO, init_dict; 
															observables=observables, kwargs...))

	function plot(P::AbstractDict)::Dict{String,Any} 
		
		plot_check_zero(P, task; 
										transf_lab="\$\\log_{10}\$",
										transf_data=apply_log10abs!
										)

	end

	return PlotTask(task, 
									[init_sliders_obs(observables);
									 [("obs_group", ChecksWLO.combs_groups()),	("ylim",)]], 
									"Curves_yofx", plot)

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

	C = Calculation("Check zeros for several $Y", ChecksWLO, init_dict; 
									observables=observables, kwargs...)

	task, out_dict, = ComputeTasks.init_multitask(C, [Y=>1], [1=>""])

	P0 = Dict("obs_group"=>"-") 

	function plot(P::AbstractDict)::Dict{String,Any} 
		
		plot_check_zero_multi(merge(P, P0), task, out_dict;
										transf_lab="\$\\log_{10}\$",
										transf_data=apply_log10abs!
										)

	end

	return PlotTask(task, 
									[init_sliders_obs(observables); ("ylim",)], 
									"Curves_yofx", 
									plot)

end 










































































































































































































































































































#############################################################################
end # module TasksPlots 
