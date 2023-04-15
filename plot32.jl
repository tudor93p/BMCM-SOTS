import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file_32.jl")

Ys = [
			:preserved_symmetries, 
			#:s0_Hamilt, :s_Hamilt, :b_Hamilt,
			]


tasks = [
				 init(BMCMSOTS, :CheckZero),  

#				 init(BMCMSOTS, :WannierBands2), 

#				 init(BMCMSOTS, :WannierBands1),

			(init(BMCMSOTS, :CheckZero_atYs; Y=Y) for Y=Ys)..., 

				init(BMCMSOTS, :CheckZero_atPS_vsX; X=:nr_kPoints), 

				init(BMCMSOTS, :WannierGap), 

				(init(BMCMSOTS,:WannierGap_atYs; Y=Y) for Y=Ys)...,

				(init(BMCMSOTS, :CheckZero_atPS_vsX_atYs; X=:nr_kPoints,Y=Y) for Y=Ys)...  ,

			];


for t in tasks[1:1]
	if !isempty(ComputeTasks.missing_data(t,show_missing=true))
		ComputeTasks.get_data_one(t,mute=false) 
	end 
end 

myPlots.plot(tasks)
