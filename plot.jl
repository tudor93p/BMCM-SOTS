import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file.jl")

tasks = [
	init(BMCMSOTS,:CheckZero),
	init(BMCMSOTS,:CheckZero_atY; Y=:preserved_symmetries),
	init(BMCMSOTS,:WannierGap),
	init(BMCMSOTS,:WannierGap_atY; Y=:preserved_symmetries),
	 ]


ComputeTasks.missing_data(tasks[1]) 

#ComputeTasks.get_data_one(task, mute=false) 

myPlots.plot(tasks)
