using Revise, Test 

import BMCMSOTS  

include("../input_file.jl")

Revise.retry()

for fn in [

"checks", 












]



	@info "********  $fn  ********"

	include("$fn.jl")


end 







nothing 
