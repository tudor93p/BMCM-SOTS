import myLibs: ComputeTasks 

import BMCMSOTS, myPlots 


include("input_file_1.jl")



tasks = [
				 init(BMCMSOTS, :WannierBands1),

			];


myPlots.plot(tasks)
