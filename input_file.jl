



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




input_checks  = Dict{Symbol,Any}(

  :allparams => (

		theta = [1/4]*2pi,
					
		nr_kPoints = [30],

		preserved_symmetries = ["All", "None", "Mx", "Ct", "Mx Ct"],
	
		nr_perturb_strength = [3],
		
		max_perturb_strength = [0.4],
	
		nr_perturb_instances = [2]

		kPoint_start = [-1],

  		), 


	:digits => (

			nr_kPoints = (3,0),

			preserved_symmetries = (),

#			perturb_strength = (1,2),
		
			nr_perturb_strength = (2,0),

			max_perturb_strength = (1,2),

			nr_perturb_instances = (2,0),

			kPoint_start = (1,2),

		),
  
	:observables => ["D110"],

	
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




