module TasksPlots
#############################################################################


#import LinearAlgebra,QuadGK,Optim

#import myLibs: Utils, ComputeTasks, Algebra, Parameters, Lattices
#import myLibs.Parameters:UODict 

import Combinatorics 


import myPlots

using myLibs.ComputeTasks: CompTask   

using myPlots: PlotTask 

import ..ChecksWLO, ..Helpers 

using ..Helpers.hParameters: Calculation  


#===========================================================================#
#
function CheckQuantiz(init_dict::AbstractDict;
											observables::AbstractVector{<:AbstractString},
											kwargs...)::PlotTask
#
#---------------------------------------------------------------------------#

	task = CompTask(Calculation(ChecksWLO, init_dict; 
															observables=observables, kwargs...))


	function plot(P::AbstractDict)::Dict{String,Any} 
		
		obs = P["obs"]

		data = task.get_data(P; fromPlot=true, mute=false, target=obs)


#		2 ["dir1","sector","stat"]
#		3 ["dir1", "stat"]  
#		2 ["stat"]
#		1 ["dir1"] 

# x  2 + 3 
# y  2 + 3

# Mx 2 
# Ct 4
# co 4 

		Y = data[obs] 

		legend  = data[ChecksWLO.fn_legend(obs)]	 

		axnames = collect(keys(legend)) 

		axvars = collect(values(legend))




		I = Vector{Union{Base.OneTo{Int},Int}}(undef, ndims(Y)-1) 

		fullax = map(in(split(P["obs_group"]," ")), axnames)
		partax = findall(!,fullax) 

		for (d,fa) in enumerate(fullax)

			fa && setindex!(I, axes(Y,d+1), d)

		end 

		ci = ChecksWLO.obs_CI(Y, P["obs_i"], partax .+ 1)
		
		for (a,i) in zip(partax,ci.I) 

			I[a] = i

		end 

		fixedlabel = ChecksWLO.get_label(ci, axvars[partax]) 

		fixedlabel = isempty(fixedlabel) ? obs : "$obs  ($fixedlabel)"


		nr_curves = mapreduce(length,*,I)


		chlabels = Vector{String}(undef,nr_curves)

		ys = [Vector{Float64}(undef,size(Y,1)) for i=1:nr_curves] 


		for (n,i) in enumerate(Iterators.product(I...))

			chlabels[n] = ChecksWLO.get_label(i[fullax], axvars[fullax]) 

			for j=axes(Y,1)

				ys[n][j] = abs(Y[j,i...])

			end 

#			copy!(ys[n], view(Y,:,i...)) 

		end  


#		y, lab = ChecksWLO.get_data_and_lab(data, obs, P["obs_i"])
														 
		x, xlabel = ChecksWLO.xxlabel(data)

		return Dict{String,Any}(

			"xs"=> [x for n=1:nr_curves],

			"ys" => ys, 

			"xlim" => extrema(x),

#			"ylim" => [-1.5,1.5],

			"labels" => chlabels,

			"ylabel" => fixedlabel, 

			"xlabel" => xlabel,

			)

	end  



	group_by = [join(c," ") for c=Combinatorics.powerset(["dir1","sector","stat"],1)]

	
	return PlotTask(task, 
									[ ("observables", observables), 
									 ("obs_index", 8), 
									 ("obs_group", group_by)
									 ],
									"Curves_yofx", plot)

end 







































































































































































































































































































#############################################################################
end # module TasksPlots 
