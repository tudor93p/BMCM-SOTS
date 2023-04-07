import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file_5.jl")


tasks = [

				init(BMCMSOTS, :CheckZero_atPS_vsX; X=:nr_kPoints), 

			];



@assert isempty(ComputeTasks.missing_data(init(BMCMSOTS, :CheckZero),show_missing=true))



myPlots.plot(tasks)
