import myLibs: ComputeTasks 
import Dates 

import BMCMSOTS  

t0 =Dates.now() 

include("input_file.jl")

tasks = [
				 init(BMCMSOTS,:CheckZero),
				 init(BMCMSOTS,:WannierBands2),
				 ]


ComputeTasks.missing_data.(tasks)
ComputeTasks.get_data_one.(tasks, mute=false) 

shuffle = gethostname()=="tudor-HP"



#@info "Preparations finished. Proceed to calculations? y/n" 

#if occursin("y",lowercase(readline(stdin))) 
if Dates.now() - t0 > Dates.Minute(5)

	ComputeTasks.get_data_all.(tasks, shuffle=shuffle, mute=false) 

end 

