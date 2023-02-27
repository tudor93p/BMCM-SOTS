import myLibs: ComputeTasks, Utils
import myPlots 

import LinearAlgebra

import BMCMSOTS:WLO, CalcWLO, ChecksWLO ,BBH

import BMCMSOTS: MB


P = (braiding_time = 0.25, kPoint_start = -1, 
		 nr_kPoints = 30,
		 preserved_symmetries = "All",#Ct",#"None",#"Ct",
		 nr_perturb_strength = 11, 
		 max_perturb_strength = 0.4, 
		 nr_perturb_instances = 1, #perturb_strength = 0.2,
		 )

@show P

function same(a,b;atol::Float64=1e-10)::Bool 

	all(<(atol), Utils.dist_periodic(a, b, 1))

end  

#function same(f::Function, a; kwargs...)::Function end  


function chapter(s::AbstractString) 

	n = 78 - 2 
	n1 = div(n-length(s),2)
	n2 = n-n1-length(s)

	println(string("\n",
								repeat("-",n1),
								" ",
								s,
								" ",
								repeat("-", n2),
								"\n"
								))

end 



println() 

observables=input_checks[:observables]
symmetries = input_checks[:allparams][:preserved_symmetries] 
	

MBparams = ChecksWLO.parse_MB_params(P)
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

perturb0_BBH = symms in ("None","All") ? perturb0 : BBH.get_perturb_on_mesh("R2", nk, k0)





#ChecksWLO.set_results_two!(results, nk, k0, 1, s1, 
#													 ChecksWLO.get_data(MB.get_psiH(MBparams, nk, k0), results) )
#

i = first(s2e)
ps = first(rest_strengths) 




#ps =  0

println("Perturbation strength: $ps")
println("Preserved symmetry: $symms\n")


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#include("wlos_MBvBBH_calc.jl")


#error() 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


chapter("Completeness: vanishing polarization")

@testset "MB+BBH: p(WLO all bands) = 0" begin 

	@test all(same.(p1all,0))
	@test all(same.(BBH_p1all,0))

end;println() 



chapter("WLO") 


@testset "MB+BBH: Σ p(wcc1) = 0" begin 

	@test same(WLO.wcc_stat(p1occup_x  + p1unocc_x), 0) 
	@test same(WLO.wcc_stat(p1occup_y  + p1unocc_y), 0) 

	@test same(WLO.wcc_stat(BBH_p1occup_x  + BBH_p1unocc_x),0) 
	@test same(WLO.wcc_stat(BBH_p1occup_y  + BBH_p1unocc_y),0)

#	@show WLO.wcc_stat(p1occup_y_  + p1unocc_y_)

end;println()



@testset "MB+BBH: Σ_i p_i(wcc1) = p(WLO all bands)" begin 

	@test same(p1occup_x  + p1unocc_x, p1all[1]) 
	@test same(p1occup_y  + p1unocc_y, p1all[2])  
	@test same(BBH_p1occup_x  + BBH_p1unocc_x, BBH_p1all[1]) 
	@test same(BBH_p1occup_y  + BBH_p1unocc_y, BBH_p1all[2])  


end;println()  

chapter("Edge polarization")


@testset "Σ_k p(wcc1,k) = 0" begin  
	
	broken=!iszero(ps)&&symms=="None"

	@test same(WLO.wcc_stat(p1occup_y)[1],0) broken=broken
	@test same(WLO.wcc_stat(p1occup_x)[1],0) broken=broken


	@test same(WLO.wcc_stat(BBH_p1occup_y)[1],0) broken=broken
	@test same(WLO.wcc_stat(BBH_p1occup_x)[1],0) broken=broken
end 



chapter("Nested WLO space decomposition")



	
tot_wcc2_x     = WLO.wcc_stat(sum(nus2pm_yx) + sum(eta2pm_yx))
tot_wcc2_y     = WLO.wcc_stat(sum(nus2pm_xy) + sum(eta2pm_xy) )
BBH_tot_wcc2_x = WLO.wcc_stat(sum(BBH_nus2pm_yx) + sum(BBH_eta2pm_yx) )
BBH_tot_wcc2_y = WLO.wcc_stat(sum(BBH_nus2pm_xy) + sum(BBH_eta2pm_xy))


@show tot_wcc2_x tot_wcc2_y 

@testset "MB: Σ_{i,k} p(wcc2,i,k) = 0 complete" begin 

		@test same(tot_wcc2_x[1],0) broken=!iszero(ps)&&symms!="All"
		@test same(tot_wcc2_y[1],0) broken=true 

	end;println()
@testset "MB: Σ_{i} p(wcc2,i) = 0 complete at each k" begin 

		@test same(tot_wcc2_x[2],0) broken=!iszero(ps)&&symms!="All"
		@test same(tot_wcc2_y[2],0) broken=true 

end;println()


@show BBH_tot_wcc2_x BBH_tot_wcc2_y
@testset "BBH: Σ_{i,k} p(wcc2,i,k) = 0 complete" begin 


	@test same(BBH_tot_wcc2_x[1],0)
	@test same(BBH_tot_wcc2_y[1],0)

end;println() 


@testset "BBH: Σ_{i} p(wcc2,i) = 0 complete at each k" begin 

	@test same(BBH_tot_wcc2_x[2],0) broken=!iszero(ps)&&symms!="All"
	@test same(BBH_tot_wcc2_y[2],0) broken=!iszero(ps)&&symms!="All"

end;println() 




#	@show WLO.wcc_stat(sum(eta2pm_xy)+sum(nus2pm_xy))
#	@show WLO.wcc_stat(sum(eta2pm_xy)-sum(nus2pm_xy))


#@show WLO.wcc_stat(sum(BBH_nus2pm_yx) + sum(BBH_eta2pm_yx))




#	broken=!iszero(ps) 
#	broken=!iszero(ps)




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

 
@testset "BBH: Σ p(wcc2) and p(wcc1)" begin  

	@show WLO.wcc_stat(sum(BBH_nus2pm_yx)-BBH_p1occup_x)

	@test same(sum(BBH_nus2pm_yx), BBH_p1occup_x) broken=!iszero(ps)&&symms!="All"

	@test same(sum(BBH_nus2pm_xy), BBH_p1occup_y) broken=!iszero(ps) &&symms!="All" 
	@show WLO.wcc_stat(sum(BBH_nus2pm_xy)-BBH_p1occup_y)

	@test same(sum(BBH_eta2pm_yx), BBH_p1unocc_x) broken=!iszero(ps) &&symms!="All"
	@test same(sum(BBH_eta2pm_xy), BBH_p1unocc_y) broken=!iszero(ps) &&symms!="All"

end;println()  



@testset "MB: Σ p(wcc2) and p(wcc1)" begin 

	@test same(sum(eta2pm_yx), p1unocc_x) broken=!iszero(ps)&&symms!="All"
	@test same(sum(nus2pm_yx), p1occup_x) broken=!iszero(ps)&&symms!="All"

	@show WLO.wcc_stat(sum(nus2pm_yx)-p1occup_x)



	@test same(sum(nus2pm_xy), p1occup_y) broken=true 

	@show WLO.wcc_stat(sum(nus2pm_xy)-p1occup_y)  

	@test same(sum(eta2pm_xy), p1unocc_y) broken=true

end;println()


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

chapter("Symmetry protection")


@testset "MB: Chiral symm" begin 

	@show WLO.wcc_stat(sum(nus2pm_xy))
	@test same(WLO.wcc_stat(sum(nus2pm_xy))[1],0) broken=true  

	@test same(WLO.wcc_stat(p1occup_y)[1], 0) 

end ;println()

@testset "MB: Mirror symm" begin 

	@test same(first.(WLO.wcc_stat.(nus2pm_xy))...)

	@test same(nus2pm_yx...)


end; println()


chapter("Nested quantization")


@testset "BBH: nested quantization x=0" begin 

	@show WLO.wcc_stat(BBH_nus2pm_yx[1])

	for BBH_nu in BBH_nus2pm_yx
		
		for BBH_mu_or_sigma in WLO.wcc_stat(BBH_nu)

			@test BBH_mu_or_sigma≈0 atol=1e-8 broken=!iszero(ps)&&symms!="All"
		
		end 

	end  

end;println()
@testset "MB: nested quantization x=0" begin 
	@show WLO.wcc_stat(nus2pm_yx[1])

	for nu in nus2pm_yx
	
		for mu_or_sigma in WLO.wcc_stat(nu) 

			@test mu_or_sigma≈0 atol=1e-7 broken=!iszero(ps)&&symms!="All"

		end 
				
	end  

end;println()

@testset "BBH: nested quantization y=0.5" begin 

	@show WLO.wcc_stat(BBH_nus2pm_xy[1],[0,0.5])   

	for BBH_nu in BBH_nus2pm_xy 

		@test WLO.wcc_stat(BBH_nu,[0,0.5])[2]≈0 atol=1e-8 broken=!iszero(ps)&&symms!="All"
		@test WLO.wcc_stat(BBH_nu,[0,0.5])[1]≈0.5 atol=1e-8 broken=!iszero(ps)&&symms!="All"


	end   

end; println() 
@testset "MB: nested quantization y=0.5" begin 

	@show WLO.wcc_stat(nus2pm_xy[1],[0,0.5])  

	for nu in nus2pm_xy 
		
		@test same(WLO.wcc_stat(nu,[0,0.5])[2],0) broken=true 
		@test same(WLO.wcc_stat(nu,[0,0.5])[1],0.5) broken=true 


	end  


end; println() 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




error() 
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
