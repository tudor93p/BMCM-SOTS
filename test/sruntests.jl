include("../input_file_32.jl") 

import BMCMSOTS 


for fn in [

#"checks", 

#"periodic_funs",

#"symms",

#"overlap_matmul",
#"wlos",

"adaptive",
#"ribbon",
#

#"wlos_MBvBBH",







]



	@info "********  $fn  ********"

	include("$fn.jl")


end 



































nothing
