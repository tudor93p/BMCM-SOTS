import myLibs: ComputeTasks, Utils
import myPlots 

import BMCMSOTS:WLO ,MB , CalcWLO, ChecksWLO 


p = (braiding_time = 0.25, nr_kPoints = 40, kPoint_start = -1, preserved_symmetries = "None", nr_perturb_strength = 6, max_perturb_strength = 0.8, nr_perturb_instances = 2, perturb_strength = 0.0)

observables=input_checks[:observables]
symmetries = input_checks[:allparams][:preserved_symmetries] 


for symm in symmetries 
	
	p1 = (braiding_time = 0.25, nr_kPoints = 40, kPoint_start = -1, preserved_symmetries = symm, nr_perturb_strength = 6, max_perturb_strength = 0.8, nr_perturb_instances = 2, perturb_strength = 0.1)

#@time data1 = ChecksWLO.Compute(p; observables=observables)

@time data2 = CalcWLO.Compute(p1; observables=observables)

end 

error() 

tasks = [
#				 init(BMCMSOTS, :CheckZero), 
				 
				 init(BMCMSOTS, :WannierBands2),
				 init(BMCMSOTS, :WannierBands1),
				 ];

pdata =map(tasks ) do task0

p = ComputeTasks.get_first_paramcomb(task0)

target = rand(observables)

P = task0.get_plotparams(p)
	P["obs_group"]= "dir1"

P["obs"] = target 
P["obs_i"] = rand(1:10)

task0.plot(P)

end 
