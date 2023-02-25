import myLibs: ComputeTasks, Utils
import myPlots 

import BMCMSOTS:WLO ,MB , CalcWLO, ChecksWLO 


P = (braiding_time = 0.25, nr_kPoints = 30, kPoint_start = -1, preserved_symmetries = "None", nr_perturb_strength = 8, max_perturb_strength = 0.6, nr_perturb_instances = 1, perturb_strength = 0.2)

observables=input_checks[:observables]
symmetries = input_checks[:allparams][:preserved_symmetries] 
	

MBtime = ChecksWLO.braiding_time(P)
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

ChecksWLO.set_results_two!(results, nk, k0, 1, s1, 
													 ChecksWLO.get_data(MB.get_psiH(MBtime, nk, k0), results) )


i = first(s2e)
ps = first(rest_strengths)

psiH = MB.get_psiH(MBtime, nk, k0, perturb0, ps)

data = ChecksWLO.get_data(psiH, results)

((eigW1_occup_x,nus2pm_xy), 
(eigW1_unocc_x,eta2pm_xy),
(eigW1_occup_y,nus2pm_yx), 
(eigW1_unocc_y,eta2pm_yx) ) = data 
		
#data_x[1] = eig1: (wf_plus, nu_plus, wf_minus, nu_minus)
#data_x[2] = nus2: (nu2_plus, nu2_minus)
#

for eig in (eigW1_occup_x, eigW1_unocc_x,
						eigW1_occup_y, eigW1_unocc_y)

	@show WLO.WannierGap_fromSubspaces(eig)

end 

p1occup_x = WLO.polariz_fromSubspaces(eigW1_occup_x)
p1unocc_x = WLO.polariz_fromSubspaces(eigW1_unocc_x)
p1occup_y = WLO.polariz_fromSubspaces(eigW1_occup_y)
p1unocc_y = WLO.polariz_fromSubspaces(eigW1_unocc_y)


p1occup_y_ = sum(nus2pm_xy) 
p1occup_x_ = sum(nus2pm_yx) 
p1unocc_y_ = sum(eta2pm_xy) 
p1unocc_x_ = sum(eta2pm_yx) 


println(LinearAlgebra.norm(Utils.dist_periodic(p1occup_y_, p1occup_y, 1)))



  





























error()

for symm in symmetries 
	
	p1 = (braiding_time = 0.25, nr_kPoints = 40, kPoint_start = -1, preserved_symmetries = symm, nr_perturb_strength = 6, max_perturb_strength = 0.8, nr_perturb_instances = 2, perturb_strength = 0.1)

#@time data1 = ChecksWLO.Compute(p; observables=observables)

@time data2 = CalcWLO.Compute(p1; observables=observables)

end 

error() 

tasks = [
#				 init(BMCMSOTS, :CheckZero), 
				 
				 init(BMCMSOTS, :WannierBands2),
				 init(BMCMSOTS, :WannierBands1),
				 ];

pdata =map(tasks ) do task0

p = ComputeTasks.get_first_paramcomb(task0)

target = rand(observables)

P = task0.get_plotparams(p)
	P["obs_group"]= "dir1"

P["obs"] = target 
P["obs_i"] = rand(1:10)

task0.plot(P)

end 
