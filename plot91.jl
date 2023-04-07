import myLibs: ComputeTasks 
import BMCMSOTS, myPlots 



include("input_file_9.jl")

B = init(BMCMSOTS)


Ys = [
			:preserved_symmetries, 
			:braiding_time,
			]


tasks = [
				B(:RibbonWannierBands1),
#				B( :RibbonWannierDensity), 
B( :RibbonWannierDensityCenters), 
B(:RibbonPolarization),
B(:RibbonPolarization_vsX; X=:braiding_time),

			];


#@assert isempty(ComputeTasks.missing_data(B(:CheckZero), show_missing=true))


for i = [1]

	if ComputeTasks.missing_data(tasks[i])>0
	
		ComputeTasks.get_data_one(tasks[i], mute=false) 
	
	end  
end 



myPlots.plot(tasks)




