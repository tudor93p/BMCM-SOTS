import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file.jl")

Ys = [:preserved_symmetries, :s0_Hamilt, :s_Hamilt, :b_Hamilt]


tasks = [
#				 init(BMCMSOTS, :CheckZero),  

				 init(BMCMSOTS, :WannierBands2),
#				 init(BMCMSOTS, :WannierBands1),
				init(BMCMSOTS, :CheckZero_atPS_vsX; X=:nr_kPoints), 

				(init(BMCMSOTS, :CheckZero_atPS_vsX_atYs; X=:nr_kPoints,Y=Y) for Y=Ys)...  ];

#tasks = vcat(
#				[init(BMCMSOTS,:CheckZero_atYs; Y=Y) for Y=Ys],
#				[init(BMCMSOTS,:WannierGap_atYs; Y=Y) for Y=Ys],
#				)


#for t in tasks 
#	@assert ComputeTasks.missing_data(t)==0
#end 

#ComputeTasks.get_data_one(tasks[1], mute=false) 

myPlots.plot(tasks)
