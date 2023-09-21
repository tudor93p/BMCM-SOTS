import BMCMSOTS 
import BMCMSOTS: FiniteSyst, MB

import myPlots  

include("../input_file_11.jl")

P = (
		 braiding_time = 1/8,
		 s0_Hamilt = 0.0,
		 s_Hamilt = 0,
		 b_Hamilt = 1,
		 width = 6,
		 )


@show MB.magnetic_field(P) MB.spin_singlet(P) 


println(round.(MB.magnetic_field(P)./maximum(abs,MB.magnetic_field(P)),digits=3))



#data = FiniteSyst.Compute(P)
#@show keys(data)

println("-----")




tasks = [
				 init(BMCMSOTS, :Spectrum0D),
				 #init(BMCMSOTS, :LocalOper0D),
#				 init(BMCMSOTS, :LocalOper0D_oneState),
#				 init(BMCMSOTS, :OperMZMs_vsX; X=:width),
#				init(BMCMSOTS, :Spectrum0D_vsX; X=:braiding_time),
				 ];


Data =map(tasks) do task1 

	pP = tasks[1].get_plotparams(P) 


	pP = Utils.adapt_merge(pP, "oper"=> "IPRq", "smooth"=>0.3)

	println() 
@show pP 

data = task1.get_data(pP,fromPlot=true,mute=false,force_comp=true)

data |> keys |> println 

pl = task1.plot(pP) 


pl |>   keys |> println

@show extrema(pl["y"])

return pl 



end 



#myPlots.plot(tasks)















