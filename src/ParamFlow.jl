# adapted from /media/tudor/Tudor/Work/2020_Snake-states/SnakeStates/Helpers 
#
module hParameters
#############################################################################

#using Constants: DATAROOT 

import ...DATAROOT 

import myLibs: Parameters, Utils
using myLibs.Parameters: UODict, ODict



export Calculation


#===========================================================================#
#
#
#---------------------------------------------------------------------------#

function with_deps(M::Module)::Vector{Module}

	vcat(M, (with_deps(m) for m in Utils.getprop(M, :Dependencies, []))...)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#function find_in_modules(Ms::Utils.List, desired_k::Symbol)
#
#	for M in Ms 
#
#		M isa Module || continue 
#
#		out = Utils.getprop(M, desired_k)
#
#		isnothing(out) || return out 
#
#	end 
#
#	return nothing 
#
#end 
#
#function find_in_modules(Ms, desired_k::Symbol, output_k::Symbol)
#
#	Dict(output_k=>find_in_modules(Ms, desired_k))
#
#end 



function keep_usedkeys(M::Module)::Function
	
	uk = Parameters.union_usedkeys(with_deps(M)...)

	return function ku(d::UODict)::UODict 
		
		Parameters.keep_usedkeys(uk, d)

	end 

end 

function keep_usedkeys(M::Module, d::UODict)::UODict

	keep_usedkeys(M)(d)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function kwargs_PF(Ms::Tuple{Vararg{<:Module}},
									 args...)::Vector{Pair{Symbol,Any}}

	kwargs_PF(collect(Ms), args...)

end 

function kwargs_PF(Ms::AbstractVector{<:Module}, 
									 kwargs=[])::Vector{Pair{Symbol,Any}}

	vals = Utils.zipmap(Utils.getprop([:valid_paramcomb, :adjust_paramcomb]),
											Utils.flatmap(with_deps, Ms))

	return vcat(collect(kwargs),
							[Pair(p...) for p in zip([:constraint, :adjust], vals)])

end 


function args_PF(M::Module; 
								 allparams, 
								 digits, 
								 usedkeys=Symbol[],
								 kwargs...)::Tuple

	(Parameters.union_usedkeys(with_deps(M)..., usedkeys),
	 digits, 
	 Utils.getprop(M, :f_allparams, identity)(allparams)
	 )

end  




function args_PF(M::Module, input_dict::UODict)::Tuple 
	
	args_PF(M; input_dict...)
	
end 


function args_PF(M::Module, input_dicts::AbstractVector{<:UODict})::Tuple 
	
	args_PF(M; 
					digits = getindex.(input_dicts, :digits), 
					allparams = getindex.(input_dicts, :allparams), 
					) 


#	Utils.getprop(M, :f_allparams, identity)(

end 


function args_PF(tuples::Vararg{<:Tuple})
	
	(args_PF(t...) for t in tuples)

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function ParamFlow(NrParamSets, usedkeys; allparams::UODict,
#									 												digits::ODict,
#																					kwargs...
#									)::Parameters.ParamFlow
#									 
#	Parameters.ParamFlow(NrParamSets, allparams, usedkeys, digits, DATAROOT)
#
#end 

#function prepare_PFargs()


function ParamFlow(NrParamSets::Int, 
									 usedkeys::Union{<:Function, <:AbstractVector{<:Symbol}},
									 input_dict::UODict;
									 kwargs...
									)::Parameters.ParamFlow

	Parameters.ParamFlow(DATAROOT,
											 NrParamSets, 
											 usedkeys, 
											 input_dict[:digits], 
											 input_dict[:allparams];
											 kwargs...
											 )

end 


function ParamFlow(usedkeys::Union{<:Function, <:AbstractVector{<:Symbol}},
									 input_dict::UODict;
									 kwargs...)::Parameters.ParamFlow 

	ParamFlow(1, usedkeys, input_dict; kwargs...)

end 



function ParamFlow(M::Module, 
									 input_dict::Union{UODict, AbstractVector{<:UODict}};
									 kwargs...)::Parameters.ParamFlow

	Parameters.ParamFlow((DATAROOT,M), args_PF(M, input_dict)...; 
											 kwargs_PF([M], kwargs)...)

end 

function ParamFlow(M::Module; 
									 allparams, digits, 
									 usedkeys = [],
									kwargs...)::Parameters.ParamFlow

	Parameters.ParamFlow((DATAROOT, M),
											 args_PF(M; 
															 allparams=allparams, 
															 digits=digits,
															 usedkeys=usedkeys,
															 );
											 kwargs_PF([M], kwargs)...)

end




function ParamFlow(tuples::Vararg{<:Tuple}; kwargs...)

	i = findlast(t->isa(t[1], Module), tuples)

	isnothing(i) && return Parameters.ParamFlow(DATAROOT, 
																							args_PF(tuples...)...; 
																							kwargs...)
		

	return Parameters.ParamFlow((DATAROOT, tuples[i][1]), 
															args_PF(tuples...)...;
															kwargs_PF(first.(tuples[1:i]), kwargs)...)

end 

#path, tuple_level1, tuple_level2, ...
#(DATAROOT,M),tuples...





#path, NrParamSets, tuple_level1, tuple_level2, ...






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


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function merge_input_dicts(dicts::Vararg{<:UODict})::Dict 

	merge(merge, Dict(), dicts...)

end 

function merge_input_dicts(list_dicts::Utils.List, dicts::Vararg{<:UODict})::Vector{<:UODict}

	[merge_input_dicts(d, dicts...) for d in list_dicts]

end 




function merge_allparams_1(lead_params::UODict; kwargs...)

	merge_input_dicts(kwargs, Dict(:allparams=>lead_params), ) 

end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Calculation(M::Module, input_dict::UODict; kwargs...)

	Parameters.Calculation(ParamFlow(M, input_dict), M; kwargs...)

end 


function Calculation(M::Module; kwargs...)

	Parameters.Calculation(ParamFlow(M; kwargs...), M; kwargs...)

end 

function Calculation(name::AbstractString, M::Module; kwargs...)

	Parameters.Calculation(name, 
												 ParamFlow(M; kwargs...), M; kwargs...)

end 

function Calculation(name::AbstractString,
										 M::Module, input_dict::UODict; kwargs...)

	Parameters.Calculation(name,
												 ParamFlow(M, input_dict), M; kwargs...)

end 


function Calculation(name::AbstractString,
										 tuples::Vararg{Tuple}; kwargs...)

	i = findlast(t->isa(t[1], Module), tuples)::Int
	
	return Parameters.Calculation(name,
																ParamFlow(tuples...; kwargs...), 
																tuples[i][1]; kwargs...)

end 

function Calculation(tuples::Vararg{Tuple}; kwargs...)

	i = findlast(t->isa(t[1], Module), tuples)::Int
	
	return Parameters.Calculation(ParamFlow(tuples...; kwargs...), 
																tuples[i][1]; kwargs...)

end 






























#############################################################################
end # module hParameters 


