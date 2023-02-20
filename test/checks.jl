import Combinatorics

import myLibs: ComputeTasks, Utils
import myPlots 

task = init(BMCMSOTS,:CheckZero)


#ComputeTasks.missing_data(task);

#ComputeTasks.get_data_one(task, mute=false); 
#ComputeTasks.get_data_one(task, Utils.DictRandVals; mute=false);  


ComputeTasks.get_data_all(task; check_data=false, mute=false); 

#for P in task.get_paramcombs()
##
##	target = [obs for obs=input_checks[:observables] if !task.files_exist(P, target=obs)]
#
#	target= rand(input_checks[:observables])
#
#	for target in input_checks[:observables]
#
##	@show target 
#
#	y = task.get_data(P; target=target, mute=true)[target]
#
#
#@assert 	BMCMSOTS.ChecksWLO.all_symms_preserved(P...)==all(selectdim(y,1,1)≈selectdim(y,1,i) for i=axes(y,1))
#
#end 
#	
#
#end









#for g in ["-";[join(c," ") for c=Combinatorics.powerset(["dir1","sector","stat"],1)]]

#task_y = init(BMCMSOTS,:CheckZero_atY; Y=:nr_kPoints) 

#for obs in input_checks[:observables]
#
#
#g = "dir1 sector stat"
#
#	P = task.get_plotparams(ComputeTasks.get_rand_paramcomb(task))
##	P["obs"] = rand(input_checks[:observables]) 
#
#	P["obs"]=obs 
#	P["obs_i"] = rand(1:10)
#	P["obs_group"]= g 
#
#	d = task.plot(P)
#
#	d["ylim"]
#
#	D = task_y.plot(P)
#
#
#	@show size(D["ys"])
#
#end 
#



#myPlots.plot(task,task_y)








































