import myLibs: ComputeTasks, Utils
import myPlots 

import LinearAlgebra,Statistics

import BMCMSOTS:WLO, CalcWLOadapt, ChecksWLO ,BBH


import BMCMSOTS: MB



P = (braiding_time = 0.25, 
		 s0_Hamilt = 0.1,
		 s_Hamilt = 1,
		 b_Hamilt = 1,
		 kPoint_start = -1, 
		 nr_kPoints = 20,
		 preserved_symmetries = "All",
#		 nr_perturb_strength = 3, 
#		 max_perturb_strength = 0.4, 
#		 nr_perturb_instances = 1, 
		 perturb_strength = 0.2,
		 )




CalcWLOadapt.Compute(P; observables=input_checks[:observables])


