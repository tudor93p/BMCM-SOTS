include("../input_file.jl") 

import BMCMSOTS 


for fn in [

"checks", 

#"periodic_funs",










]



	@info "********  $fn  ********"

	include("$fn.jl")


end 



































nothing
