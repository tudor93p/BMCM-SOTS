import myLibs: ComputeTasks, Utils
import myPlots 

import LinearAlgebra

import BMCMSOTS:WLO, CalcWLO, ChecksWLO ,BBH

#import BMCMSOTS: MB


P = (braiding_time = 0.25, nr_kPoints = 60, kPoint_start = -1, preserved_symmetries = "None", nr_perturb_strength = 10, max_perturb_strength = 0.6, nr_perturb_instances = 1, perturb_strength = 0.2)

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





#psiH = MB.get_psiH(MBtime, nk, k0, perturb0, ps) 

psi_BBH = BBH.get_psiH(BBHtheta, nk, k0, perturb0, ps) 

#data = ChecksWLO.get_data(psiH, results) 

data_BBH = ChecksWLO.get_data(psi_BBH, results)


#(
#(eigW1_occup_x,nus2pm_xy), 
#(eigW1_unocc_x,eta2pm_xy),
#(eigW1_occup_y,nus2pm_yx), 
#(eigW1_unocc_y,eta2pm_yx) ) = data 
(
(BBH_eigW1_occup_x,BBH_nus2pm_xy), 
(BBH_eigW1_unocc_x,BBH_eta2pm_xy),
(BBH_eigW1_occup_y,BBH_nus2pm_yx), 
(BBH_eigW1_unocc_y,BBH_eta2pm_yx) ) = data_BBH
		
#data_x[1] = eig1: (wf_plus, nu_plus, wf_minus, nu_minus)
#data_x[2] = nus2: (nu2_plus, nu2_minus)
#

@testset "p(MB)==p(BBH)" begin 

	for (
#			 (nu2_plus,nu2_minus),
			 (BBH_nu2_plus,BBH_nu2_minus),
			 expected_nu_BBH) in zip(
#															 (nus2pm_yx, nus2pm_xy),
															 (BBH_nus2pm_yx, BBH_nus2pm_xy),
														 BBH.BBHoutcomes(BBHtheta))

#		@test Utils.closest_periodic_b(rand(nu2_minus),[0,0.5],1)≈expected_nu_BBH
		@test Utils.closest_periodic_b(rand(BBH_nu2_minus),[0,0.5],1)≈expected_nu_BBH


		#@test same(rand(BBH_nu2_minus), expected_nu_BBH) 


	end 

end 


#gaps = [WLO.WannierGap_fromSubspaces(eig) for eig in (eigW1_occup_x, eigW1_unocc_x, eigW1_occup_y, eigW1_unocc_y)]
#@show gaps 
gaps_BBH = [WLO.WannierGap_fromSubspaces(eig) for eig in (BBH_eigW1_occup_x, BBH_eigW1_unocc_x, BBH_eigW1_occup_y, BBH_eigW1_unocc_y)]

@show gaps_BBH 

println()

####

#p1occup_x = WLO.polariz_fromSubspaces(eigW1_occup_x)
#p1unocc_x = WLO.polariz_fromSubspaces(eigW1_unocc_x)
#p1occup_y = WLO.polariz_fromSubspaces(eigW1_occup_y)
#p1unocc_y = WLO.polariz_fromSubspaces(eigW1_unocc_y)

BBH_p1occup_x = WLO.polariz_fromSubspaces(BBH_eigW1_occup_x)
BBH_p1unocc_x = WLO.polariz_fromSubspaces(BBH_eigW1_unocc_x)
BBH_p1occup_y = WLO.polariz_fromSubspaces(BBH_eigW1_occup_y)
BBH_p1unocc_y = WLO.polariz_fromSubspaces(BBH_eigW1_unocc_y)

#p1occup_y_ = sum(nus2pm_xy) 
#p1occup_x_ = sum(nus2pm_yx) 
#p1unocc_y_ = sum(eta2pm_xy) 
#p1unocc_x_ = sum(eta2pm_yx) 

BBH_p1occup_y_ = sum(BBH_nus2pm_xy) 
BBH_p1occup_x_ = sum(BBH_nus2pm_yx) 
BBH_p1unocc_y_ = sum(BBH_eta2pm_xy) 
BBH_p1unocc_x_ = sum(BBH_eta2pm_yx) 

#p1xs = (p1occup_x, 
#				p1unocc_x, 
#				p1occup_x_, 
#				p1unocc_x_ )
#p1ys = (p1occup_y, 
#				p1unocc_y, 
#				p1occup_y_, 
#				p1unocc_y_ )

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

#	global p1xs,p1ys = map(enumerate((p1xs,p1ys))) do (dir1,p1s)
#
#		map(p1s) do p1 
#			
#			t,p11 = CalcWLO.check_nu_k_dep(selectdim(p1,1,1), dir1)
#
#			@test t 
#
#			return copy(p11)
#	
#		end 
#	
#	end   

	global BBH_p1xs,BBH_p1ys = map(enumerate((BBH_p1xs,BBH_p1ys))
																 ) do (dir1,p1s)

		map(p1s) do p1 
			
			t,p11 = CalcWLO.check_nu_k_dep(selectdim(p1,1,1), dir1)

			@test t 

			return copy(p11)
	
		end 
	
	end    

#	global nus2pm_yx,nus2pm_xy = map(enumerate((nus2pm_yx,nus2pm_xy))
#																	 ) do (dir2,nus2)
#
#		map(nus2) do nu2
#			
#			t,nu21 = CalcWLO.check_nu_k_dep(selectdim(nu2,1,1), dir2)
#
#			@test t 
#
#			return copy(nu21)
#	
#		end 
#	
#	end    

	global BBH_nus2pm_yx,BBH_nus2pm_xy = map(enumerate((BBH_nus2pm_yx,BBH_nus2pm_xy))) do (dir2,nus2)

		map(nus2) do nu2
			
			t,nu21 = CalcWLO.check_nu_k_dep(selectdim(nu2,1,1), dir2)

			@test t 

			return copy(nu21)
	
		end 
	
	end   

end; println() 



#p1xs,p1ys = (Tuple(copy(selectdim(selectdim(p1,1,1),d1,1)) for p1=p1s) for (d1,p1s)=enumerate((p1xs,p1ys)))

#( p1occup_x, p1unocc_x, p1occup_x_, p1unocc_x_,) = p1xs 
#( p1occup_y, p1unocc_y, p1occup_y_, p1unocc_y_,) = p1ys  

(BBH_p1occup_x, 
 BBH_p1unocc_x, 
 BBH_p1occup_x_, 
 BBH_p1unocc_x_,) = BBH_p1xs  

(BBH_p1occup_y, 
 BBH_p1unocc_y, 
 BBH_p1occup_y_, 
 BBH_p1unocc_y_,) = BBH_p1ys  


@testset "wcc1: total polariz. zero " begin 

#	@test same(p1occup_x  + p1unocc_x,0) 
#	@test same(p1occup_y  + p1unocc_y,0)
	@test same(BBH_p1occup_x  + BBH_p1unocc_x,0) 
	@test same(BBH_p1occup_y  + BBH_p1unocc_y,0)
	
	@test same(BBH_p1occup_x_  + BBH_p1unocc_x_,0) 
	@test same(BBH_p1occup_y_  + BBH_p1unocc_y_,0)
end;println() 


@testset "D127a: polariz. in two ways BBH" begin 


	@show WLO.wcc_stat(BBH_p1occup_x, [0,0.5])
	@show WLO.wcc_stat(BBH_p1unocc_x, [0,0.5])
	@show WLO.wcc_stat(BBH_p1occup_x_,[0,0.5])
	@show WLO.wcc_stat(BBH_p1unocc_x_,[0,0.5])

	println()

	@show WLO.wcc_stat(BBH_p1occup_y, [0,0.5])
	@show WLO.wcc_stat(BBH_p1unocc_y, [0,0.5])
	@show WLO.wcc_stat(BBH_p1occup_y_,[0,0.5])
	@show WLO.wcc_stat(BBH_p1unocc_y_,[0,0.5])

	println()


	diff_x_BBH = [maximum(Utils.dist_periodic(BBH_p1occup_x, p1,1)) for p1=(BBH_p1occup_x_,BBH_p1unocc_x_)]

	diff_y_BBH = [maximum(Utils.dist_periodic(BBH_p1occup_y, p1,1)) for p1=(BBH_p1occup_y_,BBH_p1unocc_y_)]

	@show diff_x_BBH diff_y_BBH
	
#	@test same(BBH_p1occup_x, BBH_p1occup_x_)|same(BBH_p1occup_x, BBH_p1unocc_x_)


end;println()


@testset "nested quantization BBH x=0" begin 

	for BBH_nu in BBH_nus2pm_yx
		
		@show WLO.wcc_stat(BBH_nu,[0,0.5])  
				
	end  

end 

@testset "nested quantization BBH y=0.5" begin 

	for BBH_nu in BBH_nus2pm_xy 
		
		@show WLO.wcc_stat(BBH_nu,[0,0.5])  

	end  


end; println() 


error() 



@testset "D127a: polariz. in two ways" begin 

	@show same(p1occup_x, p1occup_x_)
	@show same(p1occup_x, p1unocc_x_)

	diff_x = [maximum(Utils.dist_periodic(p1occup_x, p1,1)) for p1=(p1occup_x_,p1unocc_x_)]
	diff_y = [maximum(Utils.dist_periodic(p1occup_y, p1,1)) for p1=(p1occup_y_,p1unocc_y_)]

	@show diff_x diff_y 


end;println()
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
