using Revise, Test 

import BMCMSOTS  

include("../input_file.jl")



for fn in [

"checks", 












]



	@info "********  $fn  ********"

	include("$fn.jl")


end 







nothing 
