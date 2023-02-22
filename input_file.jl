import JLD, OrderedCollections 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




input_checks  = Dict{Symbol,Any}(

  :allparams => (

		braiding_time = [1/4],
					
		nr_kPoints = 10:20:210,

		kPoint_start = [-1], 

		preserved_symmetries = ["None", "Mx", "My", "Ct", "Mx+Ct", "My+Ct", "All"],
	
		nr_perturb_strength = [21],
		
		max_perturb_strength = [0.8],
	
		nr_perturb_instances = [10],

  		), 


	:digits => (

			braiding_time = (1,2),

			nr_kPoints = (3,0),

			kPoint_start = (1,2), 

			preserved_symmetries = (),

			nr_perturb_strength = (2,0),

			max_perturb_strength = (1,2),

			nr_perturb_instances = (2,0),


		),
  
	:observables => ["D110","D111","D48","D30","D123","D125","D127a","WannierGap"],

	
)





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




