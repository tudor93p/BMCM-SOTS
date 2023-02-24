import myLibs: ComputeTasks, Utils
import myPlots 

import BMCMSOTS:WLO ,MB , CalcWLO 


P = (braiding_time = 0.25, nr_kPoints = 10, kPoint_start = -1, preserved_symmetries = "None", nr_perturb_strength = 21, max_perturb_strength = 0.8, nr_perturb_instances = 10, perturb_strength = 0.2) 



data = CalcWLO.Compute(P)

