import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file.jl")

Ys = [:preserved_symmetries, :nr_kPoints, :s0_Hamilt, :s_Hamilt,:b_Hamilt]



tasks = vcat(
				[init(BMCMSOTS,:CheckZero_atY; Y=Y) for Y=Ys],
				[init(BMCMSOTS,:WannierGap_atY; Y=Y) for Y=Ys],
				)


#for t in tasks 
#	@assert ComputeTasks.missing_data(t)==0
#end 

#ComputeTasks.get_data_one(tasks[1], mute=false) 

myPlots.plot(tasks)
