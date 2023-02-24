import myLibs: ComputeTasks, Utils
import myPlots 

import BMCMSOTS:WLO ,MB , CalcWLO, ChecksWLO 


P = (braiding_time = 0.25, nr_kPoints = 10, kPoint_start = -1, preserved_symmetries = "None", nr_perturb_strength = 21, max_perturb_strength = 0.8, nr_perturb_instances = 10, perturb_strength = 0.2) 


observables=input_checks[:observables]


data1 = ChecksWLO.Compute(P; observables=observables)


data2 = CalcWLO.Compute(P; observables=observables)

