#!/home/tudor/.juliaup/bin/julia




path = [
				"/mnt/Work/2018_Higher-Order-Topology/codes/BMCMSOTS/",
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
#							 "taranis" => 1, 
#							 "neper" => 1, 
#							 "horon" => 1, 
#							 "sia" => 1, 
							 "shu" => 1, 
							 "re" => 1, 
#							 "toad" => 4,
#							 "spaceml4" =>1,
#							 "nut"=>2,
#							 "yoshi" =>2,
#							 "kis" =>4,
							),
					Dict(
#							 "kis" => 6, 
#							 "spaceml4" => 4, 
							 "nut" => 8,
							 ),
					Dict(
							 "yoshi" => 12, 
#							 "toad" => 24, 

							 ),

					Dict("toad"=>3,
							 "nut"=>1,
							 "kis"=>1,
							 "spaceml4"=>1,
							 )
					]



tiers4 = [
					Dict(
							 "toad"=>47,
							 "yoshi"=>25,#60,
							 "spaceml4"=>7,
							 "nut"=>4,
							 "kis"=>4,
							 "horon"=>1,
#							 "taranis"=>1,
							 "neper"=>1,
							 "sia"=>1,
							 "re"=>1,
							 "uneg"=>1,
)
					]



tiers = tiers4
#tiers = tiers2
#tiers=tiers3 

println("HINT:  julia launch_jobs/monitor_workstation.jl "*join(union(keys.(tiers)...)," "))



LaunchJobs.getrun_commands(ARGS, tiers, path)


#LaunchJobs.getrun_commands(ARGS, [Dict("nut"=>20)], path)




