






#===========================================================================#
#
# Calculate wavefunctions and WLOs
#
#---------------------------------------------------------------------------#




#psiH = MB.get_psiH(MBparams, nk, k0, perturb0, ps) 

@time begin 
	
	global psiH, energies = MB.get_enpsiH(MBparams, nk, k0, perturb0, ps) 

#@assert psiH≈psiH_ 

#psi_BBH = BBH.get_psiH(BBHtheta, nk, k0, perturb0, ps)  
	global psi_BBH, energies_BBH = BBH.get_enpsiH(BBHtheta, nk, k0, perturb0_BBH, ps)  
	global data = ChecksWLO.get_data(psiH, results) 

	global data_BBH = ChecksWLO.get_data(psi_BBH, results)

end 


#@assert psi_BBH≈psi_BBH_ 

@testset "Energy gaps >0.8" begin 

	@test all(>(0.8), selectdim(energies,1,3)-selectdim(energies,1,2))
	@test all(>(0.8), selectdim(energies_BBH,1,3)-selectdim(energies_BBH,1,2))

end 


(
(eigW1_occup_x,nus2pm_xy), 
(eigW1_unocc_x,eta2pm_xy),
(eigW1_occup_y,nus2pm_yx), 
(eigW1_unocc_y,eta2pm_yx) ) = data 
(
(BBH_eigW1_occup_x,BBH_nus2pm_xy), 
(BBH_eigW1_unocc_x,BBH_eta2pm_xy),
(BBH_eigW1_occup_y,BBH_nus2pm_yx), 
(BBH_eigW1_unocc_y,BBH_eta2pm_yx) ) = data_BBH 


eigW1_all = [WLO.Wannier_subspaces_on_mesh!(WLO.wlo1_on_mesh(d,psiH)) for d=1:2]

BBH_eigW1_all = [WLO.Wannier_subspaces_on_mesh!(WLO.wlo1_on_mesh(d,psi_BBH)) for d=1:2]



#data_x[1] = eig1: (wf_plus, nu_plus, wf_minus, nu_minus)
#data_x[2] = nus2: (nu2_plus, nu2_minus)
#

p1all = [sum(eachslice(WLO.polariz_fromSubspaces(e),dims=1)) for e=eigW1_all]
BBH_p1all = [sum(eachslice(WLO.polariz_fromSubspaces(e),dims=1)) for e=BBH_eigW1_all]

p1occup_x = WLO.polariz_fromSubspaces(eigW1_occup_x)
p1unocc_x = WLO.polariz_fromSubspaces(eigW1_unocc_x) 
p1occup_y = WLO.polariz_fromSubspaces(eigW1_occup_y)
p1unocc_y = WLO.polariz_fromSubspaces(eigW1_unocc_y)

BBH_p1occup_x = WLO.polariz_fromSubspaces(BBH_eigW1_occup_x)
BBH_p1unocc_x = WLO.polariz_fromSubspaces(BBH_eigW1_unocc_x)
BBH_p1occup_y = WLO.polariz_fromSubspaces(BBH_eigW1_occup_y)
BBH_p1unocc_y = WLO.polariz_fromSubspaces(BBH_eigW1_unocc_y)





#p1xs = (p1occup_x, 
#				p1unocc_x, 
#				p1occup_x_, 
#				p1unocc_x_ )
#
#p1ys = (p1occup_y, 
#				p1unocc_y, 
#				p1occup_y_, 
#				p1unocc_y_ )
#
#BBH_p1xs = (BBH_p1occup_x, 
#				BBH_p1unocc_x, 
#				BBH_p1occup_x_, 
#				BBH_p1unocc_x_ ) 
#
#BBH_p1ys = (BBH_p1occup_y, 
#				BBH_p1unocc_y, 
#				BBH_p1occup_y_, 
#				BBH_p1unocc_y_ ) 
#
#println(LinearAlgebra.norm(Utils.dist_periodic(p1occup_y_, p1occup_y, 1)))
#===========================================================================#
#
# Easy sanity checks
#
#---------------------------------------------------------------------------#

chapter("Easy sanity checks")

@testset "corner polarization around 0 and 1/2" begin 

	for (
			 (nu2_plus,nu2_minus),
			 (BBH_nu2_plus,BBH_nu2_minus),
			 expected_nu_BBH) in zip(
															 (nus2pm_yx, nus2pm_xy),
															 (BBH_nus2pm_yx, BBH_nus2pm_xy),
														 BBH.BBHoutcomes(BBHtheta))

		@test Utils.closest_periodic_b(rand(nu2_minus),[0,0.5],1)≈expected_nu_BBH
		@test Utils.closest_periodic_b(rand(BBH_nu2_minus),[0,0.5],1)≈expected_nu_BBH


		#@test same(rand(BBH_nu2_minus), expected_nu_BBH) 


	end 

end; println()


gaps = [WLO.WannierGap_fromSubspaces(eig) for eig in (eigW1_occup_x, eigW1_unocc_x, eigW1_occup_y, eigW1_unocc_y)]

gaps_BBH = [WLO.WannierGap_fromSubspaces(eig) for eig in (BBH_eigW1_occup_x, BBH_eigW1_unocc_x, BBH_eigW1_occup_y, BBH_eigW1_unocc_y)]

@testset "Wannier gaps >1e-2" begin 

	@test all(>(1e-2),gaps)
	@test all(>(1e-2),gaps_BBH)

end; println()





@testset "BBH: wcc1 and wcc2 dep. on perp. k only" begin 
	global BBH_p1s = global (
					(BBH_p1occup_x, BBH_p1unocc_x, ),
					(BBH_p1occup_y, BBH_p1unocc_y, )
					) = map(enumerate((
										(BBH_p1occup_x, BBH_p1unocc_x, ),
										(BBH_p1occup_y, BBH_p1unocc_y, )
														 ))
																 ) do (dir1,p1s)

		map(p1s) do p1 
			
			t,p11 = CalcWLO.check_nu_k_dep(selectdim(p1,1,1), dir1, onewarning=true)

			@test t 

			return copy(p11)
	
		end 
	
	end     


	global BBH_p1all = map(enumerate(BBH_p1all)) do (dir1,p1)

		t,p11 = CalcWLO.check_nu_k_dep(p1, dir1, onewarning=true)

		@test t 

		return copy(p11) 

	end 



	global (BBH_nus2pm_yx, BBH_nus2pm_xy,
					BBH_eta2pm_yx, BBH_eta2pm_xy,
					)= map(enumerate((
														BBH_nus2pm_yx, BBH_nus2pm_xy, 
														BBH_eta2pm_yx, BBH_eta2pm_xy, 
													 ))) do (dir2,nus2)


		map(nus2) do nu2
			
			t,nu21 = CalcWLO.check_nu_k_dep(selectdim(nu2,1,1), mod(dir2-1,2)+1,
																			onewarning=true)

			@test t 

			return copy(nu21)
	
		end 
	
	end   

end; println() 




@testset "MB: wcc1 and wcc2 dep. on perp. k only" begin 

	global p1s = global (
					( p1occup_x, p1unocc_x, ),
					( p1occup_y, p1unocc_y, )
					) = map(enumerate((
													( p1occup_x, p1unocc_x, ),
													( p1occup_y, p1unocc_y, )
														 ))) do (dir1,p1s)

		map(p1s) do p1 
			
			t,p11 = CalcWLO.check_nu_k_dep(selectdim(p1,1,1), dir1,
																		 onewarning=true)

			@test t 

			return copy(p11)
	
		end 
	
	end   


	global p1all = map(enumerate(p1all)) do (dir1,p1)

		t,p11 = CalcWLO.check_nu_k_dep(p1, dir1; onewarning=true)

		@test t 

		return copy(p11) 

	end 

#end;println()
#
#@testset "MB: wcc2 on perp. k only" begin 

	global (nus2pm_yx,
					nus2pm_xy,
					eta2pm_yx,
					eta2pm_xy,
					) = map(enumerate((
														 nus2pm_yx,
														 nus2pm_xy, 
														 eta2pm_yx,
														 eta2pm_xy, 
														 ))
																	 ) do (dir2,nus2)

		map(nus2) do nu2
			
			t,nu21 = CalcWLO.check_nu_k_dep(selectdim(nu2,1,1), mod(dir2-1,2)+1,
																			onewarning=true)

			@test t 

			return copy(nu21)
	
		end 
	
	end    

end;println() 






#p1occup_y_ = sum(nus2pm_xy) 
#p1occup_x_ = sum(nus2pm_yx) 
#p1unocc_y_ = sum(eta2pm_xy) 
#p1unocc_x_ = sum(eta2pm_yx) 

#BBH_p1occup_y_ = sum(BBH_nus2pm_xy) 
#BBH_p1occup_x_ = sum(BBH_nus2pm_yx) 
#BBH_p1unocc_y_ = sum(BBH_eta2pm_xy) 
#BBH_p1unocc_x_ = sum(BBH_eta2pm_yx) 




