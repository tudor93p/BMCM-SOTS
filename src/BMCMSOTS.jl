module BMCMSOTS
#############################################################################

const MBPATH = if gethostname()=="tudor-HP" 
	
		"/media/tudor/Tudor/Work/2018_Higher-Order-Topology/BMCMSOTS"

									else 

		"/net/horon/scratch/pahomit/BMCM-SOTS"

									end  

const PATH_SNAKE = if gethostname()=="tudor-HP" 
	
				"/media/tudor/Tudor/Work/2020_Snake-states/SnakeStates"

									else 

				"/net/horon/scratch/pahomit/SnakeStates"

									end 


const DATAROOT = "$MBPATH/Data" 

const FILE_STORE_METHOD = "jld"







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

include("Helpers.jl")

include("WLO.jl") 

include("MB.jl")
include("BBH.jl")



include("CalcWLO.jl")

include("ChecksWLO.jl")

include("TasksPlots.jl")






































































































#############################################################################
end # module BMCMSOTS
