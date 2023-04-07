import myLibs: ComputeTasks 
import BMCMSOTS, myPlots 



include("input_file_9.jl")

B = init(BMCMSOTS)


Ys = [
			:preserved_symmetries, 
			:braiding_time,
			]


tasks = [
 B( :CheckZero),  
B(:CheckZero_atYs; Y=:preserved_symmetries), 
B( :CheckZero_atPS_vsX; X=:nr_kPoints  ), 
B( :CheckZero_atPS_vsX; X=:braiding_time),
B( :CheckZero_atPS_vsX_atYs; X=:braiding_time, Y=:nr_kPoints),
B( :CheckZero_atPS_vsX_atYs; Y=:braiding_time, X=:nr_kPoints),
			];





for i = [1]

	if ComputeTasks.missing_data(tasks[i])>0

		ComputeTasks.get_data_one(tasks[i], mute=false) 
	
	end   

end 
nothing 


myPlots.plot(tasks)




