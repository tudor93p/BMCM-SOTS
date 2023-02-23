import myLibs: ComputeTasks 

import BMCMSOTS  


include("input_file.jl")

task = init(BMCMSOTS,:CheckZero) 


ComputeTasks.missing_data(task)
ComputeTasks.get_data_one(task, mute=false) 



println("Preparations finished. Proceed to calculations? y/n") 


if occursin("y",lowercase(readline(stdin))) 

	ComputeTasks.get_data_all(task, mute=false) 

end 

