import Combinatorics, LinearAlgebra ,Statistics 


import myLibs: ComputeTasks, Utils
import myPlots 
 
import BMCMSOTS:WLO ,MB , ChecksWLO

import BMCMSOTS.Helpers:Symmetries 

@testset "halfspace indices" begin 
	symm_list = filter(input_checks[:allparams][:preserved_symmetries]) do S

		for s in ["All","None","+"]

			occursin(s,S) && return false 

		end 

		return true 

	end 
	@show symm_list 

	L = 12
	for p =1:L, q=1:L
	
		for n in [2p,2p+1], m in [2q,2q+1]
		
			k0= m *pi/(n-1) 
		
			i1 = WLO.uniqueInds_kMirror(n,k0) 
			
			i2 = WLO.ind_minusk.(i1,n,k0) 
	
			@test Set(union(i1,i2))==Set(1:n-1) 
			
			for i in intersect(i1,i2) 
				@test WLO.ind_minusk(i,n,k0)==i 
			end 


			for s in symm_list 

				fij = MB.getOp_fij(s,n,k0)
				ui = Set(Base.product(MB.getOp_uniqueInds(s,n,k0)...)) 

				union!(ui, [fij(I) for I in ui])

				ti = Set(Base.product(1:n-1,1:n-1))


				@test isempty(setdiff(ui,ti))
				@test isempty(setdiff(ti,ui))
				if !isempty(setdiff(ti,ui))

					@show length(ti) length(ui)

					@show n k0 s
				@show Set(Base.product(MB.getOp_uniqueInds(s,n,k0)...)) 

				end 



			end # s in symm_list 
		

			end #nm  
			
	
	end # pq 
end # testset 



@testset "Perturb. on mesh" begin 

	symm_list = filter(!in(("All","None")),input_checks[:allparams][:preserved_symmetries])

	for n=5:13, k0=rand(0:10,3).*pi/(n-1) 

		psi0 = ChecksWLO.get_psiH(pi/2, n, k0) 

		for symms in symm_list 

#			@show symms 

			pert = ChecksWLO.get_perturb_on_mesh(symms, n, k0, 3)  

			for p in pert 
				
				@test all(ChecksWLO.has_symm_on_mesh(p, symms, n, k0))
#				if !all(ChecksWLO.has_symm_on_mesh(p, symms, n, k0))
#			
#					@show LinearAlgebra.norm(p) 
#
#					@show symms 
#
#
#
#					error() 
#
#					c = ChecksWLO.has_symm_on_mesh(p, symms, n, k0) 
#
#					for i=1:n-1,j=1:n-1 
#
#						c[1,i,j] && continue 
#
#						println((i,j))
#
#						fij = MB.getOp_fij(symms, n, k0)
#						kij = WLO.get_kij(n,k0)
#
#						println(fij(i,j))
#
#						pij1 = WLO.select_mesh_point(p,i,j)  
#						pij2 = WLO.select_mesh_point(p,fij(i,j))
#
#
#						ksum = LinearAlgebra.norm(Utils.bring_periodic_to_interval.(kij(i,j) + kij(fij(i,j)),-pi,pi))
#
#						@test ksum<1e-12 
#
#
#						@show LinearAlgebra.norm(MB.operTt(pij1)-pij2 )
#
#
#						println() 
#
#
#					end 
#
#
#				
#					@show symms 
#
#					error()
#				end  
#
#
				for symms2 in symm_list 

					occursin("+",symms2) && continue 
					occursin("+",symms) && continue 

					isempty(setdiff(MB.sepSymmString(symms2),
													MB.sepSymmString(symms))) && continue 

					@test !any(ChecksWLO.has_symm_on_mesh(p, symms2, n, k0))

				end 

				
				for s in [0,1e-8, 1e-3]

					psi1 = ChecksWLO.get_psiH(pi/2, n, k0, p, s)
	
					function diff01(arg) 
	
						ov = abs.(WLO.overlap(WLO.select_mesh_point([psi1,psi0],arg...)...)) 
	
						ov .-= one(ov) 
	
						return LinearAlgebra.norm(ov)
	
					end 
					
					d = Statistics.mean(diff01, Base.product(1:n-1,1:n-1)) 
	
	#				@show (s,d)
		end 
		end 

		end 
	end 

end 


P = (braiding_time = 0.25, nr_kPoints = 10, kPoint_start = -1, preserved_symmetries = "None", nr_perturb_strength = 21, max_perturb_strength = 0.8, nr_perturb_instances = 10)


data = ChecksWLO.Compute(P; observables=input_checks[:observables])





task0 = init(BMCMSOTS,:CheckZero)




#tasks = [
#				task0,
#				 init(BMCMSOTS,:CheckZero_atY; Y=:preserved_symmetries),
#				 init(BMCMSOTS,:WannierGap),
#				 init(BMCMSOTS,:WannierGap_atY; Y=:preserved_symmetries),
#				 ]


ComputeTasks.missing_data(task0);

#ComputeTasks.get_data_one(task0, mute=false); 
#ComputeTasks.get_data_one(task0, Utils.DictRandVals; mute=false);  


#ComputeTasks.get_data_all(task0; check_data=false, mute=false); 

p = ComputeTasks.get_first_paramcomb(task0)

#for p in task0.get_paramcombs()[1:1]

	target = rand(input_checks[:observables])

#	for target in input_checks[:observables]

#	@show target 

#	y = task0.get_data(P; target=target, mute=true)[target] 

	P = task0.get_plotparams(p)

	P["obs"] = target 
	P["obs_i"] = rand(1:10)
	P["obs_group"]= "dir1"
	P["obs_group"] = rand(BMCMSOTS.ChecksWLO.combs_groups())

	task0.plot(P)

#	for t in tasks 
#
#		t.plot(P)
#
#	end 
	

#end





#for g in ["-";[join(c," ") for c=Combinatorics.powerset(["dir1","sector","stat"],1)]]


#for obs in input_checks[:observables]
#
#
#g = "dir1 sector stat"
#
#
#	d = task0.plot(P)
#
#	d["ylim"]
#
#	D = task_y.plot(P)
#
#
#	@show size(D["ys"])
#
#end 
#



#myPlots.plot(tasks)








































