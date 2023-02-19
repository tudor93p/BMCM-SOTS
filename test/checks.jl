import Combinatorics

import myLibs: ComputeTasks, Utils
import myPlots 

task = init(BMCMSOTS,:CheckQuantiz)


#ComputeTasks.missing_data(task);

#ComputeTasks.get_data_one(task, mute=false); 
#ComputeTasks.get_data_one(task, Utils.DictRandVals; mute=false);  


#ComputeTasks.get_data_all(task; check_data=true, mute=false); 

#for P in task.get_paramcombs()
#
#	target = [obs for obs=input_checks[:observables] if !task.files_exist(P, target=obs)]
#
#	@show target 
#
#	task.get_data(P; target=target, mute=false)
#
#end
#








for g in ["-";[join(c," ") for c=Combinatorics.powerset(["dir1","sector","stat"],1)]]

	P = task.get_plotparams(ComputeTasks.get_rand_paramcomb(task))
	P["obs"] = rand(input_checks[:observables])
	P["obs_i"] = rand(1:10)
	P["obs_group"]= g 

	task.plot(P)

end 


#
#myPlots.plot(task)
#
#
#
#
#



































