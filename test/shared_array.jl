using Distributed 

@everywhere using Distributed 
@everywhere using DistributedArrays, SharedArrays, LinearAlgebra
import BMCMSOTS: WLO , ChecksWLO 

#import Random 

include("mock_H.jl")



psi1 = WLO.psiH_on_mesh(50, 0 , H) 
@show typeof(psi1)
psi2 = WLO.psiH_on_mesh(50, 0 ,H; parallel=true)  
@show typeof(psi2)

@testset "simualate ind distrib" begin 
	
	
	A = WLO.init_storage(H(rand(2)), 50; parallel=true)
	B = WLO.inds_distrib_array(H(rand(2)), 50; parallel=true)
	@show A.indices[:] B 

	@test Set(A.pids)==Set(keys(B))
@show A.pids  

for (i,p) in zip(A.indices,A.pids)

		@test i==B[p] 
	end 

	map(zip(A.indices,A.pids)) do (i,p) 

		@spawnat p begin 

			@assert DistributedArrays.localindices(A)==B[p]==i 
			println("tested $p")

		end   

	end  .|> fetch  




end 

psi3 = WLO.psiH_on_mesh(50, 0 ,H; parallel=true, shared=true)

@show typeof(psi3)

#
#pert = (rand(ComplexF64,size(psi1)),) 
#
#psi1 = WLO.psiH_on_mesh(50, 0, pert..., H) 
#psi2 = WLO.psiH_on_mesh(50, 0, pert..., H; parallel=true) 
#psi3 = WLO.psiH_on_mesh(50, 0, pert..., H; parallel=true, shared=true)
#

#println()
#
#@time "psi no distr" psi1 = WLO.psiH_on_mesh(500, 0 , H) 
#@time "psi distr" psi2 = WLO.psiH_on_mesh(500, 0 ,H; parallel=true)
#@time "psi shared" psi3 = WLO.psiH_on_mesh(500, 0 ,H; parallel=true,shared=true)
#
#println()
#


@testset "psi same" begin 

	@show typeof(psi1) typeof(psi2) typeof(psi3)
#	@test psi2 isa Array  
#	@test psi1 isa Array 
#	@show norm(psi1)
#	@show norm(psi2)

	@test norm(psi1) > 1e-10 
	@test norm(psi1-psi2) <1e-10 
	@test norm(psi1-psi3) <1e-10 


end 

println() 



out_single = WLO.get_wlo_data_mesh(psi1, true, 2, true)
@show norm.(out_single)

@show typeof.(out_single)

println()  

psi2 = convert(Array, psi2)

out_multi = WLO.get_wlo_data_mesh(psi2, true, 2, true; parallel=true)
@show norm.(out_multi)
@show typeof.(out_multi)

println()  

out_shared = WLO.get_wlo_data_mesh(psi3, true, 2, true; parallel=true,
																	 shared=true)
@show norm.(out_shared)
@show typeof.(out_shared)

@time "w1+w2 single" WLO.get_wlo_data_mesh(psi1, true, 2, true)
@time "w1+w2 multi"  WLO.get_wlo_data_mesh(psi2, true, 2, true; parallel=true)
@time "w1+w2 shared"  WLO.get_wlo_data_mesh(psi3, true, 2, true; parallel=true,shared=true)
#
@testset "wlo same" begin 

#	@test w1_multi isa Array  
#	@test w1_single isa Array 

for (S,M,SS) in zip(out_single,out_multi,out_shared)

		for (s,m,ss) in zip(S,M,SS) 
	
			@test norm(s) > 1e-10 
			@test norm(ss) > 1e-10 
			@test norm(m) > 1e-10 
			@test norm(s-m) <1e-10  
			@test norm(s-ss) <1e-10  
	
		end 
end 

end 



P = (braiding_time = 0.25, 
		 s0_Hamilt = 0.1, 
		 s_Hamilt = 1.0, 
		 b_Hamilt = 1, 
		 nr_kPoints = 20,#00, 
		 kMesh_model = "Uniform",
		 kPoint_start = -1, 
		 preserved_symmetries = "All", 
		 nr_perturb_strength = 11, 
		 max_perturb_strength = 0.6, 
		 nr_perturb_instances = 1, 
		 perturb_strength = 0.3,
		 ) 



obs=["D48"]

@time data1 = ChecksWLO.Compute_(P, obs; observables=obs) 


