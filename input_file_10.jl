import JLD, OrderedCollections, Combinatorics 

import myLibs: Utils 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

interdep_symms = Set{String}[Set(["Tt", "P", "Ct"]), Set(["Mx", "TC2y", "Tt"]), Set(["Mx", "TC2y", "P", "Ct"])] 

function good_comb(c::AbstractVector{<:AbstractString})::Bool 

	for s in interdep_symms

		issubset(s,c) && return false 

	end 

	return true 

end 
obs_unocc = ["D110", "D111", "D123", "D125"]



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




input_checks  = Dict{Symbol,Any}(

  :allparams => (

		braiding_time = [1/4],

		s0_Hamilt = [0.1],#0.3],
		
		s_Hamilt = [1],#1.5,0.5],

		b_Hamilt = [1],
					
		nr_kPoints = vcat(
											10,
											50,
											95,
											145, 
											215, 
											300,
											601,
											1001,
											2001,
											3001, # >nut 
											6001,
 										10001, 
											25001,
											40_001,
																																															 ),



		kMesh_model = ["Uniform"],#"Line"],#["Square","Sin", "Expminv"],


		kPoint_start = [-1], 

		preserved_symmetries = ["All";
#														"None";
#														[join(c,"+") for c=Combinatorics.powerset(["P", "Mx", "Ct", "Tt","TC2y"],1) if good_comb(c)][13:end];
#														"None"
														],

#		preserved_symmetries = ["None", "P", "Mx", "Ct", "Tt", "All"],


		nr_perturb_strength = [11],
		
		max_perturb_strength = [0.6],

		nr_perturb_instances = [1],

		perturb_strength = [0],#0:0.2:0.8,
#		perturb_strength =  vcat(Utils.uniqlogsp(1e-4,0.2,8,3; Trunc=true), 0.3:0.1:0.6),


  		), 


	:digits => (

			braiding_time = (1,2),

			s0_Hamilt = (1,2), 

			s_Hamilt = (1,2),

			b_Hamilt = (1,2),



			nr_kPoints = (3,0),

			kPoint_start = (1,2), 

			kMesh_model = (),



			preserved_symmetries = (),

			nr_perturb_strength = (2,0),

			max_perturb_strength = (1,2),

			nr_perturb_instances = (2,0),
			
			perturb_strength = (1,3),


		),
  
	:observables => setdiff!([
										"D110","D111","D113",
									 "D48","D30","D123",
									 "D125","D127","D127.1",
									 "WannierGap",
									 "WannierBands1",
									 "WannierBands2",
									 ],
					obs_unocc)

	
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

	println()
	
	@info "********  $f  ********"

	include("$f.jl")

end 




println("Input file loaded\n")


nothing 




