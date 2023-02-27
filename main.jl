import myLibs: ComputeTasks 

import BMCMSOTS  


include("input_file.jl")

tasks = [
				 init(BMCMSOTS,:CheckZero),
				 init(BMCMSOTS,:WannierBands2),
				 ]


ComputeTasks.missing_data.(tasks)
ComputeTasks.get_data_one.(tasks, mute=false) 



@info "Preparations finished. Proceed to calculations? y/n" 

if occursin("y",lowercase(readline(stdin))) 

	ComputeTasks.get_data_all.(tasks, shuffle=true, mute=false) 

end 

