include("../input_file.jl") 

import BMCMSOTS 


for fn in [

#"checks", 

#"periodic_funs",


"wlos",







]



	@info "********  $fn  ********"

	include("$fn.jl")


end 



































nothing
