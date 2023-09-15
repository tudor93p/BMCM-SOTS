import BMCMSOTS: FiniteSyst  

import myPlots  

include("../input_file_11.jl")

P = (
		 braiding_time = 1/8,
		 s0_Hamilt = 0.0,
		 s_Hamilt = 1,
		 b_Hamilt = 1,
		 width = 5,
		 )



#data = FiniteSyst.Compute(P)
#@show keys(data)



task1 = init(BMCMSOTS, :Spectrum0D)


pP = task1.get_plotparams(P) 

@show pP 

data = task1.get_data(pP,fromPlot=true,mute=false) 

data |> keys |> println 







myPlots.plot(task1)















