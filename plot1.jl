import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


#include("input_file_11.jl")
include("input_file_14.jl")



tasks = [
				 init(BMCMSOTS, :Spectrum0D),
				 init(BMCMSOTS, :OperMZMs_vsX; X=:width),
			];


myPlots.plot(tasks)
