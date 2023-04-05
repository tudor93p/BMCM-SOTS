import myLibs: ComputeTasks, Utils
import myPlots 

import LinearAlgebra,Statistics

import BMCMSOTS:WLO, CalcWLO, ChecksWLO ,BBH

import BMCMSOTS: MB


P = (braiding_time = 0.25, 
		 s0_Hamilt = 0.1,
		 s_Hamilt = 1,
		 b_Hamilt = 1,
		 kPoint_start = -1, 
		 nr_kPoints = 10,
		 preserved_symmetries = "All",#Ct",#"None",#"Ct",
		 nr_perturb_strength = 3, 
		 max_perturb_strength = 0.4, 
		 nr_perturb_instances = 1, 
		 perturb_strength = 0.2,
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
meshsizes = input_checks[:allparams][:nr_kPoints]
	
#
#MBparams = ChecksWLO.parse_MB_params(P)
#BBHtheta = -pi/4 
#
#nk = 	ChecksWLO.nr_kPoints(P)
#
#k0 = 	ChecksWLO.kPoint_start(P)
#
#symms = ChecksWLO.preserved_symmetries(P)
#
#strengths = ChecksWLO.perturb_strengths(P)
#
#zero_strength, rest_strengths = Base.Iterators.peel(strengths)
#s1, s2e = Base.Iterators.peel(eachindex(strengths))
#
#@assert iszero(zero_strength)
#
#
#results = ChecksWLO.init_results(strengths, ChecksWLO.get_target(observables; observables=observables))
#
#trials = ChecksWLO.get_perturb_on_mesh(P, 3268)  
#
#
#j = 1 
#
#perturb0 = trials[j] 
#
#perturb0_BBH = symms in ("None","All") ? perturb0 : BBH.get_perturb_on_mesh("R2", nk, k0)
#
#
#
#
#i = first(s2e)
#ps = first(rest_strengths) 




##ps =  0
#
#println("Perturbation strength: $ps")
#println("Preserved symmetry: $symms\n")
#

#for symm in symmetries 
#	
#	p1 = (braiding_time = 0.25, nr_kPoints = 40, kPoint_start = -1, preserved_symmetries = symm, nr_perturb_strength = 6, max_perturb_strength = 0.8, nr_perturb_instances = 2, perturb_strength = 0.1)
#
#
#@time data2 = CalcWLO.Compute(P; observables=observables)


for n in meshsizes 

	n>60 && break 

	println() 

	@show n 

	p1 = Utils.adapt_merge(P, :nr_kPoints=>n)

	data1 = ChecksWLO.Compute(p1; observables=observables) 

	println() 

	break 

end 


#
#end 
#
#
tasks = [
				 init(BMCMSOTS, :CheckZero), 
				init(BMCMSOTS, :CheckZero_atYs; Y=:preserved_symmetries),
#				 init(BMCMSOTS, :WannierBands2),
##				 init(BMCMSOTS, :WannierBands1),
#				init(BMCMSOTS, :CheckZero_atPS_vsX; X=:nr_kPoints),
#				init(BMCMSOTS, :CheckZero_atPS_vsX_atYs; X=:nr_kPoints,Y=:s_Hamilt),
#init(BMCMSOTS, :CheckZero_atPS_vsX_atYs; X=:nr_kPoints,Y=:preserved_symmetries),
#				init(BMCMSOTS, :CheckZero_atPS_vsX_atYs; X=:nr_kPoints,Y=:b_Hamilt),
#				init(BMCMSOTS, :CheckZero_atPS_vsX_atYs; X=:nr_kPoints,Y=:s0_Hamilt),
				 ];

task0 = init(BMCMSOTS, :CheckZero);

pdata =map(tasks ) do task

p = ComputeTasks.get_first_paramcomb(task0)

#@assert p[1]==task0.get_paramcombs()[1][1]

#@show p[1]

target = rand(observables)
target = "D48"

P = task0.get_plotparams(p)

#P["perturb_strength"] = 0.4

P["obs_group"]= "dir1"
P["obs_group"]= "sector"

P["obs"] = target 
P["obs_i"] = rand(1:10)
P["obs_i"] = 1 

d = task.plot(P)

@show Statistics.mean(Statistics.mean,d["ys"])
@show d["ylabel"] d["labels"]
return d 

end  



