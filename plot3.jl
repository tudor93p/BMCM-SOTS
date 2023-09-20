import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file_13.jl")



tasks = [
				 init(BMCMSOTS, :Spectrum0D),
				 init(BMCMSOTS, :Spectrum0D_vsX; X=:braiding_time),
			];


myPlots.plot(tasks)
