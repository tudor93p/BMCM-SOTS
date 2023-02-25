include("../input_file.jl") 

import BMCMSOTS 


for fn in [

#"checks", 

#"periodic_funs",

#"symms",

"wlos",







]



	@info "********  $fn  ********"

	include("$fn.jl")


end 



































nothing
