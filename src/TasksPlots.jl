module TasksPlots
#############################################################################


#import LinearAlgebra,QuadGK,Optim

#import myLibs: Utils, ComputeTasks, Algebra, Parameters, Lattices
#import myLibs.Parameters:UODict 




import myLibs: Utils, ComputeTasks
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

function apply_log10abs(src::AbstractVector{<:Real}
												)::Vector{Float64}

	apply_log10abs!(Vector{Float64}(undef, length(src)), src)

end 




function apply_log10abs!(dst::AbstractVector{Float64},
											src::AbstractVector{<:Real}
											)::AbstractVector{Float64}

	m::Float64 = -35 

	for (j,s) in enumerate(src)

		setindex!(dst, log10(abs(s)), j)

		if isfinite(dst[j]) 
			
			dst[j]<m && @warn "Update the minimum: value $(dst[j]) encountered"

		else 

			setindex!(dst, m, j)

		end 

	end 

	return dst 

end  

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
function CheckZero(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	task = CompTask(Calculation("Check zeros",
															ChecksWLO, init_dict; 
															observables=observables, kwargs...))


	function plot(P::AbstractDict)::Dict{String,Any} 
		
		data = task.get_data(P; fromPlot=true, mute=false, target=P["obs"])

#		y1, lab1 = ChecksWLO.extract_data_and_lab(data, obs, P["obs_i"];
#																							transf_data=apply_log10abs)


		ys, chlabels, fixedlab = ChecksWLO.extract_data_and_lab(data, P["obs"], P["obs_group"], P["obs_i"]; transf_data=apply_log10abs!)


		fixedlab = myPlots.join_label("\$\\log_{10}\$"*P["obs"],
																	(isempty(fixedlab) ? () : (fixedlab,))...)


		x, xlabel = ChecksWLO.xxlabel(data)



		return Dict{String,Any}(

			"xs"=> [x for y=ys],

			"ys" => ys, 

			"xlim" => extrema(x),

			"ylim" => get_ylim(P, ys),

			"labels" => chlabels,

			"ylabel" => fixedlab, 

			"xlabel" => xlabel,

			)

	end  

	
	return PlotTask(task, 
									[ ("observables", observables), 
									 ("obs_index", 8), 
									 ("obs_group", ChecksWLO.combs_groups()),
									 ("ylim",)
									 ],
									"Curves_yofx", plot)

end 



#===========================================================================#
#
function CheckZero_atY(init_dict::AbstractDict;
											 Y::Symbol,
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	C = Calculation("Check zeros for several $Y",
									ChecksWLO, init_dict; observables=observables, kwargs...)

	task, out_dict, = ComputeTasks.init_multitask(C, [Y=>1], [1=>""])


	@show out_dict 



	function plot(P::AbstractDict)::Dict{String,Any} 


		function apply_rightaway(data::AbstractDict, good_P
														 )
#		y1, lab1 = ChecksWLO.extract_data_and_lab(data, obs, P["obs_i"];
#																							transf_data=apply_log10abs)
			(ChecksWLO.xxlabel(data),
			 ChecksWLO.extract_data_and_lab(data, P["obs"], P["obs_group"], P["obs_i"]; transf_data=apply_log10abs!)
			 )

		end 


		data = task.get_data(P; fromPlot=true, target=P["obs"], apply_rightaway=apply_rightaway)
		
		
		curves = vcat((d[2][1] for d in data)...)

		chlabs = vcat((myPlots.join_label.(d[2][2], [string(y)]) for (d,y)=zip(data,out_dict["y"]))...)

		x,xlabel = data[1][1] 

		fixedlab = myPlots.join_label("\$\\log_{10}\$"*P["obs"],
							(isempty(data[1][2][3]) ? () : (data[1][2][3],))...,
											 out_dict["ylabel"])




		return Dict{String,Any}(

			"xs"=> [x for y=curves],

			"ys" => curves,

			"xlim" => extrema(x),

			"ylim" => get_ylim(P, curves),

			"labels" => chlabs,

			"ylabel" => fixedlab, 

			"xlabel" => xlabel,

			)

	end  

	
	return PlotTask(task, 
									[ ("observables", observables), 
									 ("obs_index", 8), 
									 ("obs_group", ChecksWLO.combs_groups()),
									 ("ylim",)
									 ],
									"Curves_yofx", plot)

end 










































































































































































































































































































#############################################################################
end # module TasksPlots 
