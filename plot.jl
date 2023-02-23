import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file.jl")

tasks = [
	init(BMCMSOTS,:CheckZero),
	init(BMCMSOTS,:CheckZero_atY; Y=:preserved_symmetries),
	init(BMCMSOTS,:CheckZero_atY; Y=:nr_kPoints),
	init(BMCMSOTS,:WannierGap),
	init(BMCMSOTS,:WannierGap_atY; Y=:preserved_symmetries),
	init(BMCMSOTS,:WannierGap_atY; Y=:nr_kPoints),
	 ]


ComputeTasks.missing_data(tasks[1]) 

ComputeTasks.get_data_one(tasks[1], mute=false) 

myPlots.plot(tasks)
