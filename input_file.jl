import JLD, OrderedCollections, Combinatorics



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




input_checks  = Dict{Symbol,Any}(

  :allparams => (

		braiding_time = [1/4],
					
		nr_kPoints = 10:20:90,

		kPoint_start = [-1], 

		preserved_symmetries = ["None"; join.(Combinatorics.powerset(["Mx", "Ct", "TC2y", "Tt"],1,2),"+"); "All"],  

		nr_perturb_strength = [25],
		
		max_perturb_strength = [0.8],

		nr_perturb_instances = [8],

		perturb_strength = 0:0.2:0.8,

  		), 


	:digits => (

			braiding_time = (1,2),

			nr_kPoints = (3,0),

			kPoint_start = (1,2), 

			preserved_symmetries = (),

			nr_perturb_strength = (2,0),

			max_perturb_strength = (1,2),

			nr_perturb_instances = (2,0),
			
			perturb_strength = (1,2),


		),
  
	:observables => ["D110","D111","D48","D30","D123","D125","D127a",
									 "WannierGap","WannierBands1","WannierBands2"],

	
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




