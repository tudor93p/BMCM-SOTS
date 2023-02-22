import Combinatorics

import myLibs: ComputeTasks, Utils
import myPlots 

task0 = init(BMCMSOTS,:CheckZero)




tasks = [
				task0,
				 init(BMCMSOTS,:CheckZero_atY; Y=:preserved_symmetries),
				 init(BMCMSOTS,:WannierGap),
				 init(BMCMSOTS,:WannierGap_atY; Y=:preserved_symmetries),
				 ]


ComputeTasks.missing_data(task0);

#ComputeTasks.get_data_one(task0, mute=false); 
#ComputeTasks.get_data_one(task0, Utils.DictRandVals; mute=false);  


#ComputeTasks.get_data_all(task0; check_data=false, mute=false); 

p = ComputeTasks.get_rand_paramcomb(task0)

#for p in task0.get_paramcombs()[1:1]

	target = rand(input_checks[:observables])

#	for target in input_checks[:observables]

#	@show target 

#	y = task0.get_data(P; target=target, mute=true)[target] 

	P = task0.get_plotparams(p)

	P["obs"] = target 
	P["obs_i"] = rand(1:10)
	P["obs_group"]= "dir1"
	P["obs_group"] = rand(BMCMSOTS.ChecksWLO.combs_groups())

	for t in tasks 

		t.plot(P)


	end 
	

#end





#for g in ["-";[join(c," ") for c=Combinatorics.powerset(["dir1","sector","stat"],1)]]


#for obs in input_checks[:observables]
#
#
#g = "dir1 sector stat"
#
#
#	d = task0.plot(P)
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



#myPlots.plot(tasks)








































