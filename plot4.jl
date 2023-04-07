import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file_4.jl")

Ys = [
			:preserved_symmetries, 
			#:s0_Hamilt, 
			:s_Hamilt, :b_Hamilt,
			]


tasks = [
				 init(BMCMSOTS, :CheckZero),  

#				 init(BMCMSOTS, :WannierBands2),
				 init(BMCMSOTS, :WannierBands1),
				init(BMCMSOTS, :CheckZero_atPS_vsX; X=:nr_kPoints), 

				(init(BMCMSOTS, :CheckZero_atPS_vsX_atYs; X=:nr_kPoints,Y=Y) for Y=Ys)...  
			];

#tasks = vcat(
#				[init(BMCMSOTS,:CheckZero_atYs; Y=Y) for Y=Ys],
#				[init(BMCMSOTS,:WannierGap_atYs; Y=Y) for Y=Ys],
#				)


for t in tasks[1:1]
#	ComputeTasks.missing_data(t)
	@assert isempty(ComputeTasks.missing_data(t,show_missing=true))

end 

for t in tasks[1:1]
#	ComputeTasks.get_data_one(t, mute=false) 
end 

myPlots.plot(tasks)
