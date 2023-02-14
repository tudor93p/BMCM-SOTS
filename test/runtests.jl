using Revise, Test 


import BMCMSOTS 

for fn in [

"checks", 












]



	@info "********  $fn  ********"

	include("$fn.jl")


end 







nothing 
