import myLibs: ComputeTasks 

import BMCMSOTS  

include("input_file.jl")

task = init(BMCMSOTS,:CheckZero) 


ComputeTasks.get_data_one(task, mute=false)


ComputeTasks.get_data_all(task, mute=false)


