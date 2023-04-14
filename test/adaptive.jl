import myLibs: ComputeTasks, Utils
import myPlots 

import LinearAlgebra,Statistics,Optim

import BMCMSOTS.WLO: select_mesh_point

import BMCMSOTS: MB

#import PyPlot 

import BMCMSOTS:WLO, CalcWLO 



#import .CalcWLO: fill_gaps_ks!, sum_kSteps_dist2pi!, sum_kSteps ,init_gaps_ks, verify_dk_bounds 

#println(CalcWLO.next_kPoint(x->-rand(), rand(), rand(),rand(2),rand()))

#error() 


#P = (braiding_time = 0.25, 
#		 s0_Hamilt = 0.1, 
#		 s_Hamilt = 1.0, #0.5, 
#		 b_Hamilt = 1, 
#		 nr_kPoints = 15, 
#		 kMesh_type = "Adaptive", 
#		 kMesh_model = "line", 
#		 kPoint_start = -1, 
#
#		 preserved_symmetries = "All", 
#		 nr_perturb_strength = 11, 
#		 max_perturb_strength = 0.6, 
#		 nr_perturb_instances = 1, 
#		 perturb_strength = 0.3,
#		 )
#
#
#
##nk = CalcWLO.nr_kPoints(P)
#k0 = CalcWLO.kPoint_start(P)
#
#
#	unif_kij = WLO.get_kij(nk, k0)
#	
#	uniq_kinds = WLO.uniqueInds_kMirror(nk,k0) 
##	ind_minusk = WLO.ind_minusk(nk,k0)
##	get_gap_at_k = (CalcWLO.calc_gap!, 
##									CalcWLO.pack_data_gap(
##													MB.get_args_psi(P),
##													(nk,k0),
##													(true,1)
##													))
##
##
##get_dk_from_gap = CalcWLO.kMesh_model(P)([0.3061018379852348, 0.5],[1.0e-7, 0.20943951023931953])
##
#results = CalcWLO.Compute(P; observables=input_checks[:observables])
#
#
#
#
#@testset "k mirror symm" begin 
#
#
#	for i in uniq_kinds 
#		j = WLO.ind_minusk(i, nk, k0)
#	
#		ki,kj = unif_kij(i,j)
#	
#		@test Utils.dist_periodic(ki,-kj,2pi)≈0 atol=1e-10
#	
#		for K in eachcol(results["ks"])
#
#			@test Utils.dist_periodic(K[i], -K[j], 2pi)≈0 atol=1e-8 #skip=true
#	
#		end
#	
#	
#	end 
#
#end 
#




task0 = init(BMCMSOTS, :CheckZero);
tasks = [
				 task0,
#				 init(BMCMSOTS, :WannierBands2),
#				 init(BMCMSOTS, :WannierBands1),
#					init(BMCMSOTS, :CheckZero_atYs; Y=:kMesh_model),
#				 init(BMCMSOTS, :CheckZero_atPS_vsX; X=:nr_kPoints),
#				 init(BMCMSOTS, :CheckZero_atPS_vsX_atYs; X=:nr_kPoints,Y=:kMesh_model),

				 ];

pdata = map(tasks) do task

p, = ComputeTasks.get_first_paramcomb(task0) 

map(["Uniform","Line"]) do mesh 

map(["All",
#		 "None",
		 ]) do s 
	println()

	@show mesh s 

	p1 = Utils.adapt_merge(p, :preserved_symmetries=>s, :kMesh_model=>mesh,
#												 :nr_kPoints=>95
												 )


	task.get_data(p1; force_comp=true, mute=true)

end 
end 

#
#p1 = task0.get_plotparams(p)
#
#p1["obs_group"]= "sector"
#p1["obs_i"] = 2 
#p1["smooth"]=0.3 
#
#
##P["kMesh_model"]="line"
#p1["kMesh_type"] = "Adaptive"
#
#d = task.plot(p1)


end 


