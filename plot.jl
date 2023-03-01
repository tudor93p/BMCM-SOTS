import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file.jl")

tasks = [
#	init(BMCMSOTS,:CheckZero_atY; Y=:preserved_symmetries),
#	init(BMCMSOTS,:CheckZero_atY; Y=:nr_kPoints),
init(BMCMSOTS,:WannierBands2),
init(BMCMSOTS,:CheckZero),
	init(BMCMSOTS,:WannierBands1),
	init(BMCMSOTS,:WannierGap),
#	init(BMCMSOTS,:WannierGap_atY; Y=:preserved_symmetries),
#	init(BMCMSOTS,:WannierGap_atY; Y=:nr_kPoints),
	 ]


for t in tasks[2:2]
#	@assert ComputeTasks.missing_data(t)==0
end 

for t in tasks[1:2]
#	ComputeTasks.get_data_one(t, mute=false) 
end 

myPlots.plot(tasks)
