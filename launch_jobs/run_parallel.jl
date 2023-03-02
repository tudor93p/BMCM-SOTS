#!/home/tudor/apps/julia/julia-1.8.5/bin/julia

Pkg.activate("..") 

import LaunchJobs

path = [
				"/media/tudor/Tudor/Work/2018_Higher-Order-Topology/BMCMSOTS/",
				"/net/horon/scratch/pahomit/BMCM-SOTS"
				]


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


tiers2 = [Dict(
							"spaceml4"=>29, 
							"yoshi"=>38, 
							"toad"=>38,
							"nut"=>23, 

							),
				 Dict(

							"re"=>11,
							"taranis"=>15,
							"horon"=>7,
							"sia"=>7,
							"neper"=>7,
							"kis"=>7,
							)
				 ]

tiers = tiers1 
tiers = tiers2

#LaunchJobs.getrun_commands(ARGS, tiers, path)


LaunchJobs.getrun_commands(ARGS, [Dict("nut"=>20)], path)




