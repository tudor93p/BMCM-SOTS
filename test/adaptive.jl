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

nk = CalcWLOadapt.nr_kPoints(P)
k0 = CalcWLOadapt.kPoint_start(P)

Hdata = MB.get_pertHdata(MB.params_fromP(P), MB.H);


WF_line = WLO.init_storage1(WLO.init_eigH(rand(2),Hdata...)[1], nk)

WF_mesh = WLO.init_storage(WLO.init_eigH(rand(2),Hdata...)[1], nk)


dir1 = 1 # WLO in x direction, spectrum dep on ky


dir2 = 3-dir1  # perp direction 



kij = WLO.get_kij(nk,k0) 

WLO.store_on_mesh!!(WLO.eigH!, nk, kij, WF_mesh, Hdata...)

wlo1_mesh = WLO.wlo1_on_mesh_inplace(dir1, WF_mesh) 

@testset "WLO on mesh or line" begin 
	
	for k2=1:nk-1 

		get_K_ = WLO.get_kij_constrained(nk, k0, 2, kij(k2)[1], dir2)

		WLO.store_on_mesh1!!(WLO.eigH!, nk, get_K_, WF_line, Hdata...)

		wlo1_line = WLO.wlo1_on_mesh1_inplace(WF_line) 
		
		for k1=1:nk-1
	
			@test WLO.select_mesh_point(wlo1_mesh, WLO.orderinds(dir1,k1,k2))≈WLO.select_mesh_point(wlo1_line, k1) 
	
	end 
	end 
	
end 

