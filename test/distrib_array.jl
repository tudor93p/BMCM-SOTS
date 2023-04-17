using Distributed 
@everywhere using DistributedArrays 
@everywhere import LinearAlgebra 
import BMCMSOTS: WLO , ChecksWLO 

import Random 


#A = dzeros(10,10)

#@distributed for i=axes(A,2)
#
#	@show localindices(A) 
#
#end 

#A = dzeros((2,2,10,10),workers(),[1,1,1,nworkers()])
#
#
#@testset "big/small inds" begin 
#	
#	inds = [sort(unique(rand(a,length(a)-3))) for a in axes(A)[3:4]] 
#
##	inds = axes(A)[3:4] 
#
#	@show inds 
#	println([q[4] for q in A.indices][:])
#	
#	
#	for (i_p,p) in enumerate(procs(A))
#	
#		result = @spawnat p begin 
#	
#			@assert myid()==p 
#			
#			i1 = localindices(A)[3]
#			i2 = localindices(A)[4]
#	
#			@show i2 #i2 
#	
#			e1 = [(i,I) for (i,I) in enumerate(indexin(i2,inds[2])) if !isnothing(I)]
#			e2 = [(i,I) for (I,i) in enumerate(indexin(inds[2],i2)) if !isnothing(i)]
#
##			e3 = filter(!isnothing∘last,enumerate(indexin(inds[2],i2)))
#
#	
#			@assert e1==e2 
#			
#			p2 = findall(!isnothing,indexin(i2,inds[2]))
#			p1 = findall(!isnothing,indexin(i1,inds[1]))
#
#			@assert all(			(p1,p2) .== findall.(!isnothing,indexin.(localindices(A)[3:4],inds)))
#
#
#
#			@show e1 p2 
#
#			for (a,b,(c,d)) in zip(p2,localindices(A)[4][p2],e1)
#
#
#				@assert a==c && b==inds[2][d]
#
#			end 
#
#
#			for (j,(J2,J))=enumerate(zip(i2,indexin(i2,inds[2]))),(i,I)=enumerate(indexin(i1,inds[1]))
#
#				isnothing(I) && continue 
#				isnothing(J) && continue 
#		
#				
#				@assert inds[1][I]== i1[i]
#				@assert inds[2][J]==i2[j]==J2 
#				
#
#				
#
#				localpart(A)[:,:,i,j] .= i1[i] + i2[j] 
#
#			
#			end  
#
#			intersect(i1,inds[1]) 
#			
##			@show indexin.(localindices(A)[end-1:end],inds)
#
#
#	
#		end  
#	
#	fetch(result)
#	println()
#	
#	
#	end 
#	
#end 


include("mock_H.jl")


@testset "findall indexin" begin 

	for n= 10:10:100

		a=sort(unique(rand(1:100,n)))
		b=sort(unique(rand(1:100,n)))
	
		@test a!=b 
	
		f1 = findall(!isnothing, indexin(a,b))
		f2 = WLO.findall_indexin(a,b)
#		f3 = WLO.findall_indexin_2(a,b)

		@test Vector(f1)==Vector(f2) # ==Vector(f3
	
	end 

end 


psi1 = WLO.psiH_on_mesh(50, 0 , H) 
psi2 = WLO.psiH_on_mesh(50, 0 ,H; parallel=true) 
pert = (rand(ComplexF64,size(psi1)),) 
#pert = ()


psi1 = WLO.psiH_on_mesh(50, 0, pert..., H) 
psi2 = WLO.psiH_on_mesh(50, 0, pert..., H; parallel=true) 



#
#@time "psi no distr" psi1 = WLO.psiH_on_mesh(100, 0 , H) 
#@time "psi distr" psi2 = WLO.psiH_on_mesh(100, 0 ,H; parallel=true)
#
@testset "psi==psi_distr" begin 

	@test psi2 isa Array  
	@test psi1 isa Array 
#	@show LinearAlgebra.norm(psi1)
#	@show LinearAlgebra.norm(psi2)

	@test LinearAlgebra.norm(psi1) > 1e-10 
	@test LinearAlgebra.norm(psi1-psi2) <1e-10 


end 

println() 

out_single = WLO.get_wlo_data_mesh(psi1, true, 2, true)
@show LinearAlgebra.norm.(out_single)

@show typeof.(out_single)

println()  

out_multi = WLO.get_wlo_data_mesh(psi2, true, 2, true; parallel=true)
@show LinearAlgebra.norm.(out_multi)
@show typeof.(out_multi)



@time "w1+w2 single" WLO.get_wlo_data_mesh(psi1, true, 2, true)
@time "w1+w2 multi"  WLO.get_wlo_data_mesh(psi2, true, 2, true; parallel=true)

@testset "wlo==wlo_distr" begin 

#	@test w1_multi isa Array  
#	@test w1_single isa Array 

	for (S,M) in zip(out_single,out_multi)

		for (s,m) in zip(S,M) 
	
			@test LinearAlgebra.norm(s) > 1e-10 
			@test LinearAlgebra.norm(s-m) <1e-10  
	
		end 
end 

end 





P = (braiding_time = 0.25, 
		 s0_Hamilt = 0.1, 
		 s_Hamilt = 1.0, 
		 b_Hamilt = 1, 
		 nr_kPoints = 500, 
		 kMesh_model = "Line", 
		 kPoint_start = -1, 
		 preserved_symmetries = "All", 
		 nr_perturb_strength = 11, 
		 max_perturb_strength = 0.6, 
		 nr_perturb_instances = 1, 
		 perturb_strength = 0.3,
		 ) 



obs=["D48"]

@time data1 = ChecksWLO.Compute_(P, obs; observables=obs) 


