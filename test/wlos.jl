import myLibs: ComputeTasks, Utils
import myPlots 

import LinearAlgebra

import BMCMSOTS:WLO, CalcWLO, ChecksWLO ,BBH

import BMCMSOTS: MB


P = (braiding_time = 0.25, kPoint_start = -1, 
		 nr_kPoints = 30, 
		 preserved_symmetries = "All", 
		 nr_perturb_strength = 10, max_perturb_strength = 0.6, nr_perturb_instances = 1, perturb_strength = 0.2)

@show P

function same(a,b;atol::Float64=1e-10)::Bool 

	all(<(atol), Utils.dist_periodic(a, b, 1))

end 

println() 

observables=input_checks[:observables]
symmetries = input_checks[:allparams][:preserved_symmetries] 
	

MBtime = ChecksWLO.braiding_time(P)
BBHtheta = -pi/4 

nk = 	ChecksWLO.nr_kPoints(P)
k0 = 	ChecksWLO.kPoint_start(P)

symms = ChecksWLO.preserved_symmetries(P)

strengths = ChecksWLO.perturb_strengths(P)

zero_strength, rest_strengths = Base.Iterators.peel(strengths)
s1, s2e = Base.Iterators.peel(eachindex(strengths))

@assert iszero(zero_strength)


results = ChecksWLO.init_results(strengths, ChecksWLO.get_target(observables; observables=observables))

trials = ChecksWLO.get_perturb_on_mesh(P, 3268) 
j = 1 

perturb0 = trials[j]

#ChecksWLO.set_results_two!(results, nk, k0, 1, s1, 
#													 ChecksWLO.get_data(MB.get_psiH(MBtime, nk, k0), results) )
#

i = first(s2e)
ps = first(rest_strengths)

println("Perturbation strength: $ps")
println("Preserved symmetry: $symms\n")




#===========================================================================#
#
# Calculate wavefunctions and WLOs
#
#---------------------------------------------------------------------------#




psiH = MB.get_psiH(MBtime, nk, k0, perturb0, ps) 

psi_BBH = BBH.get_psiH(BBHtheta, nk, k0, perturb0, ps) 

data = ChecksWLO.get_data(psiH, results) 

data_BBH = ChecksWLO.get_data(psi_BBH, results)


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


eigW1_all_x = WLO.Wannier_subspaces_on_mesh!(WLO.wlo1_on_mesh(1,psiH))
eigW1_all_y = WLO.Wannier_subspaces_on_mesh!(WLO.wlo1_on_mesh(2,psiH))

BBH_eigW1_all_x = WLO.Wannier_subspaces_on_mesh!(WLO.wlo1_on_mesh(1,psi_BBH))
BBH_eigW1_all_y = WLO.Wannier_subspaces_on_mesh!(WLO.wlo1_on_mesh(2,psi_BBH))



#data_x[1] = eig1: (wf_plus, nu_plus, wf_minus, nu_minus)
#data_x[2] = nus2: (nu2_plus, nu2_minus)
#


#===========================================================================#
#
# Easy sanity checks
#
#---------------------------------------------------------------------------#


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

end 


gaps = [WLO.WannierGap_fromSubspaces(eig) for eig in (eigW1_occup_x, eigW1_unocc_x, eigW1_occup_y, eigW1_unocc_y)]

gaps_BBH = [WLO.WannierGap_fromSubspaces(eig) for eig in (BBH_eigW1_occup_x, BBH_eigW1_unocc_x, BBH_eigW1_occup_y, BBH_eigW1_unocc_y)]

@testset "Wannier gaps >1e-2" begin 

	@test all(>(1e-2),gaps)
	@test all(>(1e-2),gaps_BBH)

end; println()


p1all_x = sum(eachslice(WLO.polariz_fromSubspaces(eigW1_all_x),dims=1))
p1all_y = sum(eachslice(WLO.polariz_fromSubspaces(eigW1_all_y),dims=1))
BBH_p1all_x = sum(eachslice(WLO.polariz_fromSubspaces(BBH_eigW1_all_x),dims=1))
BBH_p1all_y = sum(eachslice(WLO.polariz_fromSubspaces(BBH_eigW1_all_y),dims=1))

p1occup_x = WLO.polariz_fromSubspaces(eigW1_occup_x)
p1unocc_x = WLO.polariz_fromSubspaces(eigW1_unocc_x)
p1occup_y = WLO.polariz_fromSubspaces(eigW1_occup_y)
p1unocc_y = WLO.polariz_fromSubspaces(eigW1_unocc_y)

BBH_p1occup_x = WLO.polariz_fromSubspaces(BBH_eigW1_occup_x)
BBH_p1unocc_x = WLO.polariz_fromSubspaces(BBH_eigW1_unocc_x)
BBH_p1occup_y = WLO.polariz_fromSubspaces(BBH_eigW1_occup_y)
BBH_p1unocc_y = WLO.polariz_fromSubspaces(BBH_eigW1_unocc_y)

p1occup_y_ = sum(nus2pm_xy) 
p1occup_x_ = sum(nus2pm_yx) 
p1unocc_y_ = sum(eta2pm_xy) 
p1unocc_x_ = sum(eta2pm_yx) 

BBH_p1occup_y_ = sum(BBH_nus2pm_xy) 
BBH_p1occup_x_ = sum(BBH_nus2pm_yx) 
BBH_p1unocc_y_ = sum(BBH_eta2pm_xy) 
BBH_p1unocc_x_ = sum(BBH_eta2pm_yx) 

p1xs = (p1occup_x, 
				p1unocc_x, 
				p1occup_x_, 
				p1unocc_x_ )

p1ys = (p1occup_y, 
				p1unocc_y, 
				p1occup_y_, 
				p1unocc_y_ )

BBH_p1xs = (BBH_p1occup_x, 
				BBH_p1unocc_x, 
				BBH_p1occup_x_, 
				BBH_p1unocc_x_ ) 

BBH_p1ys = (BBH_p1occup_y, 
				BBH_p1unocc_y, 
				BBH_p1occup_y_, 
				BBH_p1unocc_y_ ) 

#println(LinearAlgebra.norm(Utils.dist_periodic(p1occup_y_, p1occup_y, 1)))



@testset "nu_d dep. on perp. k_d only" begin 

	global p1s = global (
					( p1occup_x, p1unocc_x, p1occup_x_, p1unocc_x_,),
					( p1occup_y, p1unocc_y, p1occup_y_, p1unocc_y_,)
					) = map(enumerate((p1xs,p1ys))) do (dir1,p1s)

		map(p1s) do p1 
			
			t,p11 = CalcWLO.check_nu_k_dep(selectdim(p1,1,1), dir1)

			@test t 

			return copy(p11)
	
		end 
	
	end   

	global (p1all_x,p1all_y) = map(enumerate((p1all_x,p1all_y))) do (dir1,p1)

		t,p11 = CalcWLO.check_nu_k_dep(p1, dir1)

		@test t 

		return copy(p11) 

	end 



	global nus2pm_yx,nus2pm_xy = map(enumerate((nus2pm_yx,nus2pm_xy))
																	 ) do (dir2,nus2)

		map(nus2) do nu2
			
			t,nu21 = CalcWLO.check_nu_k_dep(selectdim(nu2,1,1), dir2)

			@test t 

			return copy(nu21)
	
		end 
	
	end    

end;println() 

@testset "BBH nu_d dep. on perp. k_d only" begin 
	global BBH_p1s = global (
					(BBH_p1occup_x, BBH_p1unocc_x, BBH_p1occup_x_, BBH_p1unocc_x_,),
					(BBH_p1occup_y, BBH_p1unocc_y, BBH_p1occup_y_, BBH_p1unocc_y_,)
					) = map(enumerate((BBH_p1xs,BBH_p1ys))
																 ) do (dir1,p1s)

		map(p1s) do p1 
			
			t,p11 = CalcWLO.check_nu_k_dep(selectdim(p1,1,1), dir1)

			@test t 

			return copy(p11)
	
		end 
	
	end     

	global (BBH_p1all_x,BBH_p1all_y
					) = map(enumerate((BBH_p1all_x,BBH_p1all_y))) do (dir1,p1)

		t,p11 = CalcWLO.check_nu_k_dep(p1, dir1)

		@test t 

		return copy(p11) 

	end 



	global BBH_nus2pm_yx,BBH_nus2pm_xy = map(enumerate((BBH_nus2pm_yx,BBH_nus2pm_xy))) do (dir2,nus2)

		map(nus2) do nu2
			
			t,nu21 = CalcWLO.check_nu_k_dep(selectdim(nu2,1,1), dir2)

			@test t 

			return copy(nu21)
	
		end 
	
	end   

end; println() 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


@testset "wcc1: vanishing total polariz." begin 

	@test same(p1occup_x  + p1unocc_x, 0) 
	@test same(p1occup_y  + p1unocc_y, 0) 

	@test same(BBH_p1occup_x  + BBH_p1unocc_x,0) 
	@test same(BBH_p1occup_y  + BBH_p1unocc_y,0)

#	@show WLO.wcc_stat(p1occup_y_  + p1unocc_y_)

end;println()


@testset "wcc2: vanishing total polariz." begin 

	
	@test same(p1occup_x_  + p1unocc_x_, 0) 
	@test same(p1occup_y_  + p1unocc_y_, 0) broken=true  

	@test same(BBH_p1occup_x_  + BBH_p1unocc_x_,0) 
	@test same(BBH_p1occup_y_  + BBH_p1unocc_y_,0) 

	#@show WLO.wcc_stat(BBH_p1occup_y_  + BBH_p1unocc_y_)

end;println() 


@testset "D127a BBH: polariz. in two ways" begin 

	for (d1,BBH_p1d) in enumerate(BBH_p1s)

		for BBH_item1 in BBH_p1d 

			@test same(WLO.wcc_stat(BBH_item1), 0)

		end 

		diff_d1_BBH = [maximum(Utils.dist_periodic(BBH_p1d[1],item1,1)) for item1=BBH_p1d[3:end]]

#		@show diff_d1_BBH 

		@test any(Base.Fix1(same,0), diff_d1_BBH)

	end 


end;println() 


@testset "D127a: polariz. in two ways" begin 


	for (d1,p1d) in enumerate(p1s)

		for item1 in p1d 

			@test same(WLO.wcc_stat(item1), 0) skip=(d1==2)

		end 

		diff_d1 = [maximum(Utils.dist_periodic(p1d[1],item1,1)) for item1=p1d[3:end]]

#		@show diff_d1 
		@test any(Base.Fix1(same,0), diff_d1) skip=(d1==2)

	end 

end;println()




@testset "nested quantization x=0" begin 

	for BBH_nu in BBH_nus2pm_yx
		
#		@show WLO.wcc_stat(BBH_nu,[0,0.5])  
		@test WLO.wcc_stat(BBH_nu,[0,0.5])[2]≈0 atol=BBH.delta 
				
	end  
	for nu in nus2pm_yx
		
#		@show WLO.wcc_stat(nu,[0,0.5])  

		@test same(WLO.wcc_stat(nu,[0,0.5])[2],0,atol=1e-7)

				
	end  

end;println()

@testset "nested quantization y=0.5" begin 

	for BBH_nu in BBH_nus2pm_xy 
		
#		@show WLO.wcc_stat(BBH_nu,[0,0.5])   
		
		@test WLO.wcc_stat(BBH_nu,[0,0.5])[2]≈0 atol=BBH.delta 

	end   

	for nu in nus2pm_xy 
		
#		@show WLO.wcc_stat(nu,[0,0.5]) 
		@test same(WLO.wcc_stat(nu,[0,0.5])[2],0,atol=1e-4)


	end  


end; println() 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





error() 
#
#
#@testset "D125 for wcc2: total polarization zero" begin 
#
#
##	Utils.closest_periodic_shifted_a.(a,0,1)
#
#
#
#@test same(p1occup_x_ + p1unocc_x_,0)
#@test same(p1occup_y_ + p1unocc_y_,0)
#
#
#
#
#end 
#























#
#
#
#for symm in symmetries 
#	
#	p1 = (braiding_time = 0.25, nr_kPoints = 40, kPoint_start = -1, preserved_symmetries = symm, nr_perturb_strength = 6, max_perturb_strength = 0.8, nr_perturb_instances = 2, perturb_strength = 0.1)
#
##@time data1 = ChecksWLO.Compute(p; observables=observables)
#
#@time data2 = CalcWLO.Compute(p1; observables=observables)
#
#end 
#
#error() 
#
#tasks = [
##				 init(BMCMSOTS, :CheckZero), 
#				 
#				 init(BMCMSOTS, :WannierBands2),
#				 init(BMCMSOTS, :WannierBands1),
#				 ];
#
#pdata =map(tasks ) do task0
#
#p = ComputeTasks.get_first_paramcomb(task0)
#
#target = rand(observables)
#
#P = task0.get_plotparams(p)
#	P["obs_group"]= "dir1"
#
#P["obs"] = target 
#P["obs_i"] = rand(1:10)
#
#task0.plot(P)
#
#end 
