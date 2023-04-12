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
		 kMesh_model="expminv",
		 nr_kPoints = 33,
		 preserved_symmetries = "Mx",
#		 nr_perturb_strength = 3, 
#		 max_perturb_strength = 0.4, 
#		 nr_perturb_instances = 1, 
		 perturb_strength = 0.2,
		 )





nk = CalcWLOadapt.nr_kPoints(P)
k0 = CalcWLOadapt.kPoint_start(P)


	unif_kij = WLO.get_kij(nk, k0)
	
	uniq_kinds = WLO.uniqueInds_kMirror(nk,k0) 
#	ind_minusk = WLO.ind_minusk(nk,k0)
#	get_gap_at_k = (CalcWLOadapt.calc_gap!, 
#									CalcWLOadapt.pack_data_gap(
#													MB.get_args_psi(P),
#													(nk,k0),
#													(true,1)
#													))
#
#
#get_dk_from_gap = CalcWLOadapt.kMesh_model(P)([0.3061018379852348, 0.5],[1.0e-7, 0.20943951023931953])
#
results = CalcWLOadapt.Compute(P; observables=input_checks[:observables])




@testset "k mirror symm" begin 


	for i in uniq_kinds 
		j = WLO.ind_minusk(i, nk, k0)
	
		ki,kj = unif_kij(i,j)
	
		@test Utils.dist_periodic(ki,-kj,2pi)≈0 atol=1e-10
	
		for K in eachcol(results["ks"])

			@test Utils.dist_periodic(K[i], -K[j], 2pi)≈0 atol=1e-8 #skip=true
	
		end
	
	
	end 

end 









