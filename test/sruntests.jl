include("../input_file_9.jl") 

import BMCMSOTS 


for fn in [

#"checks", 

#"periodic_funs",

#"symms",

#"overlap_matmul",
#"wlos",

"ribbon",

#"wlos_MBvBBH",







]



	@info "********  $fn  ********"

	include("$fn.jl")


end 



































nothing
