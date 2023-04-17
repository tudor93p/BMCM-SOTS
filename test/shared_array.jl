using Distributed 
@everywhere using DistributedArrays, SharedArrays, LinearAlgebra
import BMCMSOTS: WLO , ChecksWLO 

#import Random 

include("mock_H.jl")



psi1 = WLO.psiH_on_mesh(50, 0 , H) 
psi2 = WLO.psiH_on_mesh(50, 0 ,H; parallel=true)  
psi3 = WLO.psiH_on_mesh(50, 0 ,H; parallel=true, shared=true)

#pert = (rand(ComplexF64,size(psi1)),) 
#
#psi1 = WLO.psiH_on_mesh(50, 0, pert..., H) 
#psi2 = WLO.psiH_on_mesh(50, 0, pert..., H; parallel=true) 
#


#
#@time "psi no distr" psi1 = WLO.psiH_on_mesh(100, 0 , H) 
#@time "psi distr" psi2 = WLO.psiH_on_mesh(100, 0 ,H; parallel=true)



@testset "psi==psi_distr" begin 

	@show typeof(psi1) typeof(psi2) typeof(psi3)
#	@test psi2 isa Array  
#	@test psi1 isa Array 
#	@show norm(psi1)
#	@show norm(psi2)

	@test norm(psi1) > 1e-10 
	@test norm(psi1-psi2) <1e-10 


end 

println() 

#out_single = WLO.get_wlo_data_mesh(psi1, true, 2, true)
#@show norm.(out_single)
#
#@show typeof.(out_single)
#
#println()  
#
#out_multi = WLO.get_wlo_data_mesh(psi1, true, 2, true; parallel=true)
#@show norm.(out_multi)
#@show typeof.(out_multi)
#
#
#
#@time "w1+w2 single" WLO.get_wlo_data_mesh(psi1, true, 2, true)
#@time "w1+w2 multi"  WLO.get_wlo_data_mesh(psi1, true, 2, true; parallel=true)
#
#@testset "wlo==wlo_distr" begin 
#
##	@test w1_multi isa Array  
##	@test w1_single isa Array 
#
#	for (S,M) in zip(out_single,out_multi)
#
#		for (s,m) in zip(S,M) 
#	
#			@test norm(s) > 1e-10 
#			@test norm(s-m) <1e-10  
#	
#		end 
#end 
#
#end 
#
#
#
#
#
#P = (braiding_time = 0.25, 
#		 s0_Hamilt = 0.1, 
#		 s_Hamilt = 1.0, 
#		 b_Hamilt = 1, 
#		 nr_kPoints = 500, 
#		 kMesh_model = "Line", 
#		 kPoint_start = -1, 
#		 preserved_symmetries = "All", 
#		 nr_perturb_strength = 11, 
#		 max_perturb_strength = 0.6, 
#		 nr_perturb_instances = 1, 
#		 perturb_strength = 0.3,
#		 ) 
#
#
#
#obs=["D48"]
#
#@time data1 = ChecksWLO.Compute_(P, obs; observables=obs) 
#
#
