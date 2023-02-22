import myLibs: ComputeTasks 

import BMCMSOTS  

#import myPlots

include("input_file.jl")

task = init(BMCMSOTS,:CheckZero) 


ComputeTasks.missing_data(task)
#ComputeTasks.get_data_one(task, mute=false)

#ComputeTasks.get_data_all(task, mute=false)

#myPlots.plot(task)
