import Dates  

t0 = Dates.now() 

import myLibs: ComputeTasks 

import BMCMSOTS  


include("input_file.jl")

tasks = [
				 init(BMCMSOTS,:CheckZero),
				 init(BMCMSOTS,:WannierBands2),
				 ];


#ComputeTasks.missing_data.(tasks)
ComputeTasks.get_data_one.(tasks, mute=false) 

#shuffle = gethostname()=="tudor-HP"



#@info "Preparations finished. Proceed to calculations? y/n" 

#if occursin("y",lowercase(readline(stdin)))  

t1 = Dates.DateTime("2023-02-27T17:33:32.010")


while Dates.now() < t1 #< Dates.Minute(9)

	sleep(1)

end 


ComputeTasks.get_data_all.(tasks, shuffle=true, seed=2, mute=true)


