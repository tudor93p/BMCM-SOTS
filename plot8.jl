import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file_8.jl")


println("Init tasks")

tasks = [
				 init(BMCMSOTS, :CheckZero),
				 init(BMCMSOTS, :WannierBands1),

				 init(BMCMSOTS, :CheckZero_atPS_vsX; X=:nr_kPoints), 
				 init(BMCMSOTS, :CheckZero_atPS_vsX_atYs; 
							X=:nr_kPoints,Y=:kPoint_start), 
#				 init(BMCMSOTS, :WannierGap_atYs; Y=:nr_kPoints),

			];



#@assert isempty(ComputeTasks.missing_data(init(BMCMSOTS, :CheckZero),show_missing=true))

for task0 in tasks[1:2]

ComputeTasks.get_data_one(task0,mute=false)

end 

for task0 in tasks[1:2]
ComputeTasks.missing_data(task0)
end 

myPlots.plot(tasks)
