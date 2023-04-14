import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file_10.jl")

Ys = [
			:s_Hamilt,
			:s0_Hamilt,
			:kMesh_model,
			]

task0 = init(BMCMSOTS, :CheckZero)

tasks = [
#				task0, 
#				 init(BMCMSOTS, :WannierBands2), 

			 init(BMCMSOTS, :WannierBands1),

#			(init(BMCMSOTS, :CheckZero_atYs; Y=Y) for Y=Ys)..., 


#				init(BMCMSOTS, :WannierGap), 

				(init(BMCMSOTS,:WannierGap_atYs; Y=Y) for Y=Ys)...,

				init(BMCMSOTS, :CheckZero_atPS_vsX; X=:nr_kPoints), 

				(init(BMCMSOTS, :CheckZero_atPS_vsX_atYs; X=:nr_kPoints,Y=Y) for Y=Ys)...  ,

			];


for t in [tasks[1],task0]
	isempty(ComputeTasks.missing_data(t,show_missing=true))
end 
#ComputeTasks.get_data_one(tasks[1], mute=false) 

myPlots.plot(tasks)
