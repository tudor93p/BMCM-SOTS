import Dates  

t0 = Dates.now() 

import myLibs: ComputeTasks 

import BMCMSOTS  


include("input_file_32.jl")

tasks = [
				 init(BMCMSOTS,:CheckZero),
				 init(BMCMSOTS,:WannierBands1),
				 ];


ComputeTasks.missing_data.(tasks,show_missing=false);



#error() 

ComputeTasks.get_data_one.(tasks, mute=false);

#shuffle = gethostname()=="tudor-HP"



#@info "Preparations finished. Proceed to calculations? y/n" 

#if occursin("y",lowercase(readline(stdin)))  


t1 = Dates.DateTime("2023-02-27T18:11:46.870")

if Dates.now()<t1 

	println("Calculations start at: ",t1)

	while Dates.now() < t1 
	
		sleep(0.1)
	
	end 

end 


ComputeTasks.get_data_all.(tasks, 
													 shuffle=true, seed=4, 
													 mute=false,
													 check_data=false,
													 )








