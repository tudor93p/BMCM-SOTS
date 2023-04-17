include("../input_file_10.jl") 

import BMCMSOTS 


for fn in [

#"checks", 

#"periodic_funs",

#"symms",

#"overlap_matmul",
#"wlos",
#"distrib_array",
#"shared_array",
#"adaptive_cluster",
"adaptive",
#"ribbon",
#

#"wlos_MBvBBH",







]



	@info "********  $fn  ********"

	include("$fn.jl")


end 



































nothing
