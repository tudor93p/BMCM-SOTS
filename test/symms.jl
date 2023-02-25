import Combinatorics, LinearAlgebra ,Statistics 


import myLibs: ComputeTasks, Utils
import myPlots 
 
import BMCMSOTS:WLO ,MB ,CalcWLO

import BMCMSOTS.Helpers:Symmetries 





@testset "halfspace indices" begin 

	symm_list = filter(input_checks[:allparams][:preserved_symmetries]) do S

		for s in ["All","None","+"]

			occursin(s,S) && return false 

		end 

		return true 

	end 

	L = 12
	for p =1:L, q=-L:L, n in [2p,2p+1], m in [2q,2q+1]  

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
		

	end # nmpq 
end # testset 



@testset "Perturb. on mesh" begin 

	symm_list = filter(!in(("All","None")),input_checks[:allparams][:preserved_symmetries])

	for n=5:13, k0=rand(-2n:2n,3).*pi/(n-1) 

		psi0 = MB.get_psiH(1/4, n, k0) 

		for symms in symm_list 

#			@show symms 

				pert = CalcWLO.get_perturb_on_mesh(symms, n, k0)

				@test all(MB.has_symm_on_mesh(pert, symms, n, k0))
#
				for symms2 in symm_list 

					occursin("+",symms2) && continue 
#					occursin("+",symms) && continue 

					s2 = MB.sepSymmString(symms2)

					s1 = MB.sepSymmString(symms) # should be larger <


					
					S2 = MB.has_symm_on_mesh(pert, symms2, n, k0)

					if issubset(s2,s1) 
						@test all(S2)
					elseif isdisjoint(s1,s2)

						if any(S2)

#							@show s1 s2 

							for i in findall(S2[1,:,:])

								for ki in WLO.get_kij(n,k0)(i.I)/pi 
								
									# only high-symm points

									@test ki≈Int(round(ki))

								end 

							end 




#							@show Statistics.mean(S2)
#	println()
#							@test false 
						else 

						end 
					end 


				end 

#				for s in [0,1e-8, 1e-3]
#
#					psi1 = MB.get_psiH(1/4, n, k0, pert, s)
#	
#					function diff01(arg) 
#	
#						ov = abs.(WLO.overlap(WLO.select_mesh_point([psi1,psi0],arg...)...)) 
#						ov .-= one(ov) 
#	
#						return LinearAlgebra.norm(ov)
#	
#					end 
#					
#					d = Statistics.mean(diff01, Base.product(1:n-1,1:n-1)) 
#	
##				@show (s,d)
#		end 

		end 
	end 

end 


