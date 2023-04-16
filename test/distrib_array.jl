using Distributed 
@everywhere using DistributedArrays 
@everywhere import LinearAlgebra 
import BMCMSOTS: WLO 

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



P = (braiding_time = 0.25, 
		 s0_Hamilt = 0.1, 
		 s_Hamilt = 1.0, 
		 b_Hamilt = 1, 
		 nr_kPoints = 15, 
		 kMesh_model = "Line", 
		 kPoint_start = -1, 
		 preserved_symmetries = "All", 
		 nr_perturb_strength = 11, 
		 max_perturb_strength = 0.6, 
		 nr_perturb_instances = 1, 
		 perturb_strength = 0.3,
		 )

@everywhere function H(ij::NTuple{2,Int},
						kpoint::Function,#Union{Function,<:AbstractArray{<:Real,4}},
						perturb::AbstractArray{<:Number,4},
					 )::Matrix{ComplexF64}

	H(kpoint(ij), view(perturb,:,:,ij...))

end 

@everywhere function H(k::AbstractVector{<:Real}, 
											 pert::AbstractMatrix{<:Number}
											 )
	a = H(k) 
	
	a += LinearAlgebra.Hermitian(pert)  

	return a 

end 

@everywhere function H((kx,ky)::AbstractVector{<:Real})

	intra = [0.902 + 0.879im 0.391 + 0.623im 0.761 + 0.283im 0.277 + 0.26im; 0.376 + 0.355im 0.946 + 0.354im 0.95 + 0.834im 0.88 + 0.583im; 0.706 + 0.241im 0.782 + 0.785im 0.161 + 0.928im 0.023 + 0.559im; 0.495 + 0.933im 0.357 + 0.076im 0.728 + 0.825im 0.873 + 0.4im] 
	inter_x = [0.935 + 0.807im 0.821 + 0.786im 0.193 + 0.631im 0.983 + 0.415im; 0.624 + 0.878im 0.531 + 0.596im 0.262 + 0.108im 0.322 + 0.023im; 0.379 + 0.303im 0.094 + 0.479im 0.894 + 0.136im 0.206 + 0.992im; 0.071 + 0.624im 0.087 + 0.871im 0.938 + 0.819im 0.687 + 0.422im] 
	inter_y = [0.795 + 0.063im 0.014 + 0.86im 0.746 + 0.544im 0.349 + 0.185im; 0.912 + 0.427im 0.19 + 0.882im 0.032 + 0.609im 0.48 + 0.974im; 0.165 + 0.523im 0.201 + 0.424im 0.504 + 0.787im 0.538 + 0.583im; 0.16 + 0.363im 0.259 + 0.906im 0.056 + 0.317im 0.657 + 0.245im]

	LinearAlgebra.Hermitian(intra+cis(kx)*inter_x+cis(ky)*inter_y)

end 



pmap(H,eachcol(rand(2,10)))

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

out_multi = WLO.get_wlo_data_mesh(psi1, true, 2, true; parallel=true)
@show LinearAlgebra.norm.(out_multi)
@show typeof.(out_multi)



@time "w1+w2 single" WLO.get_wlo_data_mesh(psi1, true, 2, true)
@time "w1+w2 multi"  WLO.get_wlo_data_mesh(psi1, true, 2, true; parallel=true)

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





