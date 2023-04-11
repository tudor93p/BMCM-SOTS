import myLibs: ComputeTasks, Utils
#import myPlots 

import LinearAlgebra,Statistics,Optim

import BMCMSOTS.WLO: select_mesh_point

import BMCMSOTS: MB

import PyPlot 

import BMCMSOTS:WLO, CalcWLOadapt, ChecksWLO ,BBH


#import .CalcWLOadapt: fill_gaps_ks!, sum_kSteps_dist2pi!, sum_kSteps ,init_gaps_ks, verify_dk_bounds 


P = (
		 braiding_time = 0.25, 
#		 braiding_time = rand(), 
		 s0_Hamilt = 0.1,
		 s_Hamilt = 1,
		 b_Hamilt = 1,
		 kPoint_start = -1, 
		 kMesh_type="Adaptive",
		 kMesh_model="line",
		 nr_kPoints = 51,
		 preserved_symmetries = "All",
#		 nr_perturb_strength = 3, 
#		 max_perturb_strength = 0.4, 
#		 nr_perturb_instances = 1, 
		 perturb_strength = 0.2,
		 )





nk = CalcWLOadapt.nr_kPoints(P)
k0 = CalcWLOadapt.kPoint_start(P)



#results = CalcWLOadapt.Compute(P; observables=input_checks[:observables])


@testset "k mirror symm" begin 

	unif = WLO.get_kij(nk, k0)
	
	for i in WLO.uniqueInds_kMirror(nk, k0)
	
		j = WLO.ind_minusk(i, nk, k0)
	
		ki,kj = unif(i,j)
	
		@test Utils.dist_periodic(ki,-kj,2pi)≈0 atol=1e-10
	
		for K in eachcol(results["ks"])

			println(Utils.dist_periodic(K[i], -K[j], 2pi))
	
		end
	
	
	end 

end 
#@show ks + reverse(ks,dims=1)













