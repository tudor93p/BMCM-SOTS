#!/home/tudor/apps/julia/julia-1.8.5/bin/julia

import LaunchJobs

path = [
				"/media/tudor/Tudor/Work/2018_Higher-Order-Topology/BMCMSOTS/",
				"/net/horon/scratch/pahomit/BMCM-SOTS"
				]


tiers = [Dict("spaceml4"=>28, 
							#"yoshi"=>12, 
							"nut"=>23, 
							"toad"=>38,
							"taranis"=>15,
							"horon"=>7,
							"sia"=>8,
							"re"=>11,
							"neper"=>8,
							"kis"=>8,
							)
				 ]



LaunchJobs.getrun_commands(ARGS, tiers, path)




