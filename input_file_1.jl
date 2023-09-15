import JLD, OrderedCollections, Combinatorics 

import myLibs: Utils 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




input_checks  = Dict{Symbol,Any}(

  :allparams => (

		Delta0 = 0.1:0.1:1,

		pairSymm = ["chiral", "helical", "ext.swave"],

		pairSymmConfig = [1,2],

		nr_kPoints = 10:10:200,

		kPoint_start = [-1], 

  		), 


	:digits => (

			pairSymm = (),

			pairSymmConfig = (1,0),

			Delta0 = (1,2),




			nr_kPoints = (3,0),

			kPoint_start = (1,2), 


			preserved_symmetries = (),

			nr_perturb_strength = (2,0),

			max_perturb_strength = (1,2),

			nr_perturb_instances = (2,0),
			
			perturb_strength = (1,3),


		),
  
	:observables => [
									 "WannierBands1",
									 ],

	
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
	args = (input_checks,)
	kwargs = input_checks 
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




