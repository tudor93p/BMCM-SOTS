#!/home/tudor/apps/julia/julia-1.8.5/bin/julia

path = [
				"/media/tudor/Tudor/Work/2018_Higher-Order-Topology/codes/BMCMSOTS/",
				"/net/horon/scratch/pahomit/BMCM-SOTS"
				]
Pkg.activate(path[1]) 


import LaunchJobs



tiers1 = [Dict(
							"spaceml4"=>29, 
							"yoshi"=>38, 
							"toad"=>38,
							"nut"=>23, 
							"re"=>11,
							"taranis"=>15,
							"horon"=>7,
							"sia"=>7,
							"neper"=>7,
							"kis"=>7,
							)
				 ]


tiers2 = [
					Dict(
							"yoshi"=>38, 
#							),
#				 Dict(

														"toad"=>40,),
				 Dict(
							"re"=>11,
							"taranis"=>15,
							"sia"=>7,
							"neper"=>7,
				"horon"=>7,
							),
				 Dict(
									"spaceml4"=>30,

							),
				 Dict(
#							"nut"=>23, 
							"nut"=>20,
							"taranis"=>12,
#							"spaceml4"=>10,
							"kis"=>8,
#							"yoshi"=>5

				),
#				 Dict(
#							"kis"=>7,
#							)
				] 

tiers3 = [Dict(
							 "re"=>4,
							 "kis"=>4,  
							 "horon"=>4, 
"nut"=>4,
"spaceml4"=>8,
							 "taranis"=>4,
								"shu"=>4,
							"neper"=>4, 
								"yoshi"=>8,

							"sia"=>4,
							"toad"=>32,
							)

					]



tiers = tiers1 
#tiers = tiers2
tiers=tiers3 


LaunchJobs.getrun_commands(ARGS, tiers, path)


#LaunchJobs.getrun_commands(ARGS, [Dict("nut"=>20)], path)




