module TasksPlots
#############################################################################


#import LinearAlgebra,QuadGK,Optim

#import myLibs: Utils, ComputeTasks, Algebra, Parameters, Lattices
#import myLibs.Parameters:UODict 


import LsqFit 


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
													)::Vector{Tuple}
													
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
tex(labs::Union{Char,<:AbstractString}...)::String = tex(string(labs...))

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

		"xlabel" => xlabel * " " * tex("\\lambda"),
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
#
#
#---------------------------------------------------------------------------#


function Rdlabel(dir1::Int)::String 

#		"R_"*
string("yx"[dir1]) 

end 

function get_eigval_number(P::AbstractDict,
													 nmax::Int,
													 default_value::T=1
													 )::Union{T,Int} where T

	haskey(P,"k") ? min(nmax,max(1,Int(round(P["k"])))) : default_value

end  


function pos_exp_val(rho::AbstractMatrix{<:Real})::Vector{Float64}

	[sum(prod,enumerate(rhoj))/sum(rhoj) for rhoj=eachcol(rho)]

end 

function pos_exp_val(rho::AbstractMatrix{<:Real},
										 d::Int)::Tuple{Vector{Float64},String}

	pos_exp_val(rho), tex("\\langle\\; ",Rdlabel(d),"\\;\\rangle")

end 

function pos_exp_val(rho::AbstractArray{<:Real,3},
										 d::Int)::Tuple{Vector{Float64},String}

	pos_exp_val(selectdim(rho,3,d),d)

end 

#===========================================================================#
#
function RibbonWannierDensityCenters(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

#	obs = "WannierDensity"
	obs = ["WannierBands1","WannierDensity"]

	@assert all(in(observables), obs)

	task = CompTask(Calculation("Ribbon Wannier density centers",
															RibbonWLO, init_dict; 
															observables=observables, kwargs...))


	function plot(P::AbstractDict)::Dict{String,Any} 

		d = get_dir1(P) 

		data = task.get_data(P; fromPlot=true, mute=false, target=obs)

#		size(data[obs]) = (nr_at, nr_wf)

		y0,ylabel = pos_exp_val(data["WannierDensity"],d)

		x = Vector(eachindex(y0))

		xlabel = "Eigenstate index"

		xlim = Utils.extend_limits(x,0.02)

		ylim = Utils.extend_limits(axes(data["WannierDensity"],1),0.02)


		out = Dict{String,Any}(
	
					 "x" => x,
	
					 "y" => y0,
						
					 "label" => "",
	
						"xlim"=> xlim,
	
						"ylim"=> ylim,
	
						"xlabel" => xlabel,
	
						"ylabel" => ylabel,

						)


		if haskey(P,"k")
		
			n = length(y0)

			j = get_eigval_number(P, n)

			out["z"] = setindex!(zeros(n), 1, n)

			out["label"] = string(j) 
			out["zlabel"] = "Selected state"

			for k in ["x","y"]

				aux = out[k][n]
				out[k][n] = out[k][j] 
				out[k][j] = aux   # move last such that it's visible 

			end 

		else 

			out["z"] = copy(selectdim(data["WannierBands1"],2,d))

			out["zlabel"] = tex(data[RibbonWLO.fn_legend("WannierBands1")]["dir"][d])

			out["zlim"] = [-0.5,0.5]

		end 
	
			return out 
	
		end 

		return PlotTask(task, [("obs_index", 2),
													 ], 
										"Scatter", plot,)
end  


#===========================================================================#
#
function RibbonWannierDensity(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	obs = "WannierDensity"

	@assert obs in observables 

	task = CompTask(Calculation("Ribbon Wannier density",
															RibbonWLO, init_dict; 
															observables=observables, kwargs...))


	function plot(P::AbstractDict)::Dict{String,Any} 

		d = get_dir1(P) 

		data = task.get_data(P; fromPlot=true, mute=false, target=obs,)


		j = get_eigval_number(P, size(data[obs],2))
		
		ylabel = tex("\\rho \\;\\;[\\;j=$j\\,]")

		y0 = copy(data[obs][:,j,d]) # Vector 

		x = Vector(eachindex(y0))

		xlabel = tex(Rdlabel(d))

		xlim = [x[1]-1,x[end]+1]


		xlim = Utils.extend_limits(x,1e-2) 

		ylim = get_ylim(0,nothing,data[obs])


		out = Dict{String,Any}(
	
					 "x" => x,
	
					 "y" => y0,
						
					 "label" => string(j),
	
						"xlim"=> xlim,
	
						"ylim"=> ylim,
	
						"xlabel" => xlabel,
	
						"ylabel" => ylabel,

#						"yline" => 0,
	
						)
	
	
			return out 
	
		end 

		return PlotTask(task, [("obs_index", 2),("choose_k",),
													 ], 
										"Curves_yofx", plot,)
end  


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_polariz(j::Nothing,
										 rho::AbstractMatrix,
										 nu::AbstractVector,
										 )::Tuple{Vector,String}

	(RibbonWLO.polarization(rho, nu),"")

end 

function get_polariz(j::Int,
										 rho::AbstractMatrix,
										 nu::AbstractVector,
										 )::Tuple{Vector,String}

	(RibbonWLO.polarization(selectdim(rho,2,j:j),selectdim(nu,1,j:j)),
					string(j,":"))

end  


function get_polariz(P::AbstractDict,
										 rho::AbstractArray,
										 nu::AbstractArray,
										 )::Tuple

	get_polariz(get_eigval_number(P, size(rho,2), nothing), rho, nu)

end  


function get_polariz(j::Union{Nothing,Int},
										 rho::AbstractArray{<:Real},
										 nu::AbstractMatrix{<:Real},
										 )::Tuple{Vector,Vector,String}

	@assert size(rho,3)==size(nu,2)==2

	(px,lab),(py,lab) = [get_polariz(j, selectdim(rho,3,d), selectdim(nu,2,d)) for d=1:2]

	return (px,py,lab)

end  

function get_polariz(P::Union{<:AbstractDict,<:Nothing,<:Int},
										 data::AbstractDict
										 )::Tuple{Vector,Vector,String}

	get_polariz(P, data["WannierDensity"], data["WannierBands1"])

end 

function get_polariz_halves(p::AbstractVector{<:Real}
														)::Vector{Float64}

	map([1:div(length(p),2), div(length(p),2)+1:length(p)]) do slice 
	
		Utils.bring_periodic_to_interval(sum(view(p, slice)), -0.5, 0.5)

	end 

end  


function f_get_polariz_halves(
															P::Union{<:AbstractDict,<:Nothing,<:Int}=nothing
															)::Function 

	function get_polariz_halves_(data::AbstractDict,
															args... 
															)::Tuple{Vector{Float64},Vector{String}}
	
		px,py,lab = get_polariz(P,data) 

		#		left: py(Rx<middle)


		return (vcat(
								 Utils.bring_periodic_to_interval(sum(px),-0.5,0.5),
								 get_polariz_halves(px),
								 Utils.bring_periodic_to_interval(sum(py),-0.5,0.5),
								 get_polariz_halves(py),
								 ),
						lab .* tex.([
										"p_x", "p_x(R_y < L_y/2)","p_x(R_y>L_y/2)",
										"p_y", "p_y(R_x<L_x/2)", "p_y(R_x>L_x/2)",
										])
								)

	end 

end 

#===========================================================================#
#
function RibbonPolarization(init_dict::AbstractDict;
#											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	task = CompTask(Calculation("Ribbon polarization",
															RibbonWLO, init_dict; 
															kwargs...))


	function plot(P::AbstractDict)::Dict{String,Any} 

		data = task.get_data(P; fromPlot=true, mute=false, target=["WannierBands1","WannierDensity"])

		d = get_dir1(P)   

		y0,lab = get_polariz(P,
													selectdim(data["WannierDensity"],3,d),
													selectdim(data["WannierBands1"],2,d),
													)

		p = trunc(sum(y0),digits=8)

		plab = "p_"*"xy"[d]

		ylabel = tex(plab,"\\,(\\,",Rdlabel(d),"\\,)")*"  [Total "*tex(plab,"=",rstrip(nr2str(p,5),'0'))*"]"



		x = Vector(eachindex(y0))

		xlabel = tex(Rdlabel(d))

		xlim = [x[1]-1,x[end]+1]

		ylim = Utils.extend_limits([-0.5, 0.5],1e-1)



		out = Dict{String,Any}(
	
					 "xs" => [x,x,x],
	
					 "ys" => [y0 .- 1, y0, y0 .+ 1],
						
					 "labels" => lab .* ["-1","0","+1"],
	
						"xlim"=> xlim,
	
						"ylim"=> ylim,
	
						"xlabel" => xlabel,
	
						"ylabel" => ylabel,

						"yline" => 0,
	
						)
	
	
			return out 
	
		end 

	return PlotTask(task, ("obs_index", 2), "Scatter", plot,)
end 
#===========================================================================#
#
function RibbonPolarization_vsX(init_dict::AbstractDict;
																X::Symbol,
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	C = Calculation("Ribbon polarization vs. $X", RibbonWLO, init_dict; 
															kwargs...)

	task, out_dict, construct_Z = ComputeTasks.init_multitask(C, [X=>1], 
																														[1=>"Edge polarization"])


	function plot(P::AbstractDict)::Dict{String,Any} 

		data = task.get_data(P; fromPlot=true, mute=false, target=["WannierBands1","WannierDensity"], apply_rightaway=f_get_polariz_halves())

		labels = data[1][2]

		ys = construct_Z(first, data)["z"]

		for j=axes(ys,2),i=axes(ys,1)

				ys[i,j] = Utils.bring_periodic_to_interval(ys[i,j],-0.25,0.75)

		end  


		ylim = [-0.25, 0.75]


		out = Dict{String,Any}(

					 "xs" => [out_dict["y"] for i=axes(ys,1)], 

					 "ys" => ys,
						
					 "labels" => labels,
	
					 "xlim"=> extrema(out_dict["y"]),
	
						"ylim"=> ylim,
	
						"xlabel" => out_dict["ylabel"],
	
						"ylabel" => "Edge polarization",

						"yline" => 0,
	
						)
	
			ComputeTasks.add_line!(out, P, X, "x") # also in plotdict_checkzero
	
			return out 
	
		end 

	return PlotTask(task, ("obs_index", 2), "Curves_yofx", plot,)

end 
#===========================================================================#
#
function RibbonWannierBands1(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	obs = ["WannierBands1","WannierDensity"]

	@assert all(in(observables),obs)


	task = CompTask(Calculation("Wannier bands ribbon",
															RibbonWLO, init_dict; 
															observables=observables, kwargs...))


	function plot(P::AbstractDict)::Dict{String,Any} 

		d = get_dir1(P) 

		data = task.get_data(P; fromPlot=true, mute=false, target=obs)  
	
		x,xlabel = RibbonWLO.xxlabel(data)
	
		legend = data[RibbonWLO.fn_legend("WannierBands1")] 
	
		y0 = selectdim(data["WannierBands1"],2,d) #Vector

		n = length(y0)

		gap = 1.0 

		for i = 1:n-1 
			
			gap = min(gap,
								Utils.reduce_dist_periodic(min, y0[i], view(y0, i+1:n), 1)
								)

		end 

		xlim = Utils.extend_limits(x,1e-2)

		xlim[2] = max(xlim[2],xlim[1]+1)

		out = Dict{String,Any}(
	
					 "xs" => [x,x,x],
 
					 "ys" => [y0 .- 1, copy(y0), y0 .+ 1],
 					
					 "labels" => ["-1","0","+1"],
	
						"xlim"=> xlim,
	
						"ylim"=> (-0.7, 0.7),
	
						"xlabel" => xlabel,
	
						"ylabel" => tex(legend["dir"][d]) * " " * Wannier_gap_text(gap),

						"yline" => 0,

						)

		if haskey(P,"k")
			
			j = get_eigval_number(P, n)

#			z = setindex!(zeros(length(y0)), 1, j)
			z = setindex!(zeros(n), 1, n)

			out["zlabel"] = "Selected state"

			out["labels"] = string(j,":") .* out["labels"]

			out["zs"] = [z,z,z]

			for k in ["xs","ys"], i=1:3 

				aux = out[k][i][n]
				out[k][i][n] = out[k][i][j] 
				out[k][i][j] = aux   # move last such that it's visible 

			end 

#			out["xline"] = j 

		else 

			z,zlab = pos_exp_val(data["WannierDensity"],d)


			out["zs"] = [z,z,z]

			out["zlabel"] = zlab

			out["zlim"] = [1,size(data["WannierDensity"],1)]

		end 
	
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



F_ifXisN = log10∘inv 

Flab_ifXisN = "\$\\log_{10}(1/N)\$"


function processConv_ifXisN!(out_dict::AbstractDict,
														 X::Symbol)::AbstractDict




	if X==:nr_kPoints 

		out_dict["x"] = F_ifXisN.(out_dict["x"])

		out_dict["xlabel"] = Flab_ifXisN 

#		out_dict["xlim"] = [0,maximum(out_dict["x"])] 

	end 

	out_dict["xlim"] = extrema(out_dict["x"]) 

	return out_dict 

end 

function xline_ifXisN!(out::AbstractDict, X::Symbol)::AbstractDict


	if X==:nr_kPoints && haskey(out,"xline")

		out["xline"] = F_ifXisN(out["xline"])

	end 

	return out 

end 

function line(x::Real, (a,b)::AbstractVector{<:Real})::Float64

	a*x+b 

end 
function line(x::AbstractVector{<:Real}, (a,b)::AbstractVector{<:Real}
						 )::Vector{Float64}

	a*x .+ b 

end 



function fit_line_conv(xdata::AbstractVector{<:Real},
											 ydata::AbstractVector{<:Real},
											 )::Tuple{Vector{Float64},String}


	if (length(xdata)!=length(ydata)) | (length(xdata)<2)
		
		return (zeros(2),"no fit")

	end 
												 
	fit = LsqFit.curve_fit(line, xdata, ydata, rand(2))

	return (LsqFit.coef(fit),
					string("slope=",nr2str(LsqFit.coef(fit)[1])))

end 

function fit_line_conv(xdata::AbstractVector{<:Real},
											 ydata::AbstractVector{<:Real},
											 start::Int,stop::Int 
											 )::Tuple{Vector{Float64},String}

#	start>0 || return (zeros(2),"no fit")

#	i = start:length(xdata) 
@show start stop 

	return fit_line_conv(
											 view(xdata,start:stop),
											 view(ydata,start:stop),
											 )
end 


function fit_line_conv(x::AbstractVector{<:Real},
											 y::AbstractVector{<:Real},
											 smooth::Float64
											 )::Tuple{Vector{Float64},String}

	0<=smooth<=1 || return (zeros(2),"no fit")

#	i = start:length(xdata) 
	#start = Int(round(smooth*(length(x)-2)+1))
	stop  = Int(round(length(x)*(1-smooth)+2*smooth))

	return fit_line_conv(x,y,
#											 start, length(x),
												1, stop,
											 )

end 

#function ys_lab_conv(xdata::AbstractVector{<:Real},
#										 args...
#										 )::Tuple{Vector{Float64},String}
#
#	ab,lab = fit_line_conv(xdata, args...)
#
#	return line(xdata,ab), string("slope=",nr2str(ab[1]))
#
#end 


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

	processConv_ifXisN!(out_dict, X)

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

		n = length(chlabs)

		nr_curves = n* (1+(X==:nr_kPoints))

		labels = vcat((chlabs for i=0:(X==:nr_kPoints))...)


		ys = [zeros(length(data)) for i=1:nr_curves]

		for (i,d) in enumerate(data)
			for (j,y) in enumerate(d[3])

				ys[j][i]=y

			end 
		end 


#		ys = mapreduce(Base.Fix2(getindex,3),hcat,data)


		ylabel = myPlots.join_label(transf_lab*P["obs"],
																data[1][5], 
																tex("\\lambda="*nr2str(data[1][1],2)),
																						 )

		if X==:nr_kPoints 


			for j=1:n
	
				ab,lab = fit_line_conv(out_dict["x"], ys[j], 
															 Float64(get(P,"smooth",0)))

#				y,lab = ys_lab_conv(out_dict["x"], ys[j], s)

				setindex!(ys, line(out_dict["x"],ab), j+n)
				
				labels[j+n] = myPlots.join_label(labels[j+n], lab)
	
			end 

		end 





		out = Dict{String,Any}(

							"xs" => [out_dict["x"] for i=1:nr_curves],
							
							"ys" => ys,

							"xlim" => out_dict["xlim"],

							"ylim" => get_ylim(P, view(ys, 1:n)),

							"labels" => labels,

							"ylabel" => ylabel,

							"xlabel" => out_dict["xlabel"],

						)


#										transf_lab="\$\\log_{10}\$",
#										transf_data=apply_log10abs!
#										)

		ComputeTasks.add_line!(out, P, X, "x") # also in plotdict_checkzero

		return xline_ifXisN!(out, X)

	end

	return PlotTask(task,  
									Tuple[init_sliders_obs(observables_);
									 [("obs_group", ChecksWLO.combs_groups()),	
									 ("ylim",), 
									 (X==:nr_kPoints ? (("smoothen",),) : ())...,
									 ]
									 ], 
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

	processConv_ifXisN!(out_dict, X)

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
																tex("\\lambda="*nr2str(data[1][1],2)),
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

	return xline_ifXisN!(out, X)

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
