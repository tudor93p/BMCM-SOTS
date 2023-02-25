import JLD, OrderedCollections, Combinatorics 

import myLibs: Utils 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




input_checks  = Dict{Symbol,Any}(

  :allparams => (

		braiding_time = [1/4],
					
		nr_kPoints = 10:20:70,

		kPoint_start = [-1], 

		preserved_symmetries = ["None"; join.(Combinatorics.powerset(["Mx", "Ct", "Tt"],1,2),"+"); "All"],  

		nr_perturb_strength = [8],
		
		max_perturb_strength = [0.6],

		nr_perturb_instances = [1],

#		perturb_strength = 0:0.2:0.8,
		
		perturb_strength =  vcat(Utils.uniqlogsp(1e-4,0.2,8,3; Trunc=true),
														 0.4:0.2:0.6),

  		), 


	:digits => (

			braiding_time = (1,2),

			nr_kPoints = (3,0),

			kPoint_start = (1,2), 

			preserved_symmetries = (),

			nr_perturb_strength = (2,0),

			max_perturb_strength = (1,2),

			nr_perturb_instances = (2,0),
			
			perturb_strength = (1,3),


		),
  
	:observables => ["D110","D111","D48","D30","D123","D125","D127a",
									 "WannierGap",
									 "WannierBands1",
									 "WannierBands2",
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




