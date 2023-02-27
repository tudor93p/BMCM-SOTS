import myLibs: ComputeTasks, Utils
import myPlots 

import LinearAlgebra

import BMCMSOTS:WLO, CalcWLO, ChecksWLO ,BBH

import BMCMSOTS: MB


P = (braiding_time = 0.25, s0_Hamiltonian = 0.1,
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
@time data2 = CalcWLO.Compute(P; observables=observables)
@time data1 = ChecksWLO.Compute(P; observables=observables)
#
#end 
#
#error() 
#
tasks = [
				 init(BMCMSOTS, :CheckZero), 
				 
				 init(BMCMSOTS, :WannierBands2),
				 init(BMCMSOTS, :WannierBands1),
				 ];
#
pdata =map(tasks ) do task0

p = ComputeTasks.get_first_paramcomb(task0)

target = rand(observables)

P = task0.get_plotparams(p)
	P["obs_group"]= "dir1"

P["obs"] = target 
P["obs_i"] = rand(1:10)

task0.plot(P)

end  



