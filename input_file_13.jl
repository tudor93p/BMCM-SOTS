import JLD, OrderedCollections, Combinatorics 

import myLibs: Utils 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




input_dict  = Dict{Symbol,Any}(

  :allparams => (


		 braiding_time = Utils.uniqlinsp(1/8,3/8,57,3,Trunc=true),

		s0_Hamilt = [0.0],
		
		s_Hamilt = [0],
		
		b_Hamilt = [1],

		width = vcat(
								 5,
								 25,
#								 40,
								),
				
  		), 


	:digits => (

			braiding_time = (1,3),

			s0_Hamilt = (1,2), 

			s_Hamilt = (1,2),

			b_Hamilt = (1,2),


			width = (3,0),

		),

	:operators => String["IPR"],


	
)





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





function init( target_module::Module, target_function::Symbol ;
							other_kwargs...)
###
	args = (input_dict,)
	kwargs = input_dict 
###

	task = getproperty(target_module.TasksPlots, target_function)

	return task(args...; kwargs..., other_kwargs...)

end 


function init( target_module::Module; kwargs1...)::Function

	function init_(target_function::Symbol; kwargs2...)
		
		init(target_module, target_function; kwargs1..., kwargs2...)

	end 

end 




function pr_in(f::AbstractString) 
	
	@info "********  $f  ********"

	include("$f.jl")

end 







nothing 




