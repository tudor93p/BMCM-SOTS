#include("../input_file_9.jl") 

import myPlots 
import BMCMSOTS: RibbonWLO , MB , WLO

import myLibs:TBmodel,Lattices 

observables=RibbonWLO.calc_observables  

for ps = LinRange(0,1,15)

P = (width=60,
			braiding_time = 0.25, 
			s0_Hamilt = 0.1, s_Hamilt = 1, b_Hamilt = 1, 
			kPoint_start = -1, nr_kPoints = 150, 
#			preserved_symmetries = "All",#Ct",#"None",#"Ct",  
#			preserved_symmetries = "Ct",
			preserved_symmetries = "Ct",
			nr_perturb_strength = 3, max_perturb_strength = 0.4, 
			nr_perturb_instances = 1, perturb_strength = ps,) 


data = RibbonWLO.Compute(P; observables=observables)["WannierBands1"]

nu = [data[1,1], data[end,1]]

println(ps,"\t",findmin(Utils.dist_periodic(nu,0.5,1)))

end 






#	nk = RibbonWLO.nr_kPoints(P) 
#
#	w = RibbonWLO.width(P) 
#
#	k0 = RibbonWLO.kPoint_start(P)
#
#	results = RibbonWLO.init_results(w, observables)
#
#
#	data = map(1:2) do dir1
#
#		latt = RibbonWLO.lattice(w, 3-dir1)
#	
#		h = Matrix∘TBmodel.Bloch_Hamilt(Lattices.NearbyUCs(latt); nr_orb=4, Hopping=MB.get_hoppf(P))
#
#		psi = WLO.psiH_on_mesh1(nk, k0, h)
#
#		wbb, nu = WLO.get_wlo_data_mesh1(psi) 
#
#		return psi,wbb,nu
#
##		set_results_onedir!(results, dir1, data)
#
#end  
#


#
#tasks = [
# init(BMCMSOTS, :RibbonWannierBands1),
# init(BMCMSOTS, :RibbonWannierDensityCenters),
# init(BMCMSOTS, :RibbonWannierDensity),
#init(BMCMSOTS, :RibbonPolarization) 
#]
#
##myPlots.plot(t2, t4, t3, t) 
#
#pdata = map(tasks) do t 
#
#map([merge(t.get_plotparams(P),Dict("k"=>2.3)),
#		 delete!(t.get_plotparams(P),"k")]) do P1 
#
#t.plot(P1)
#
#	end 
#end 
#
#t0 = init(BMCMSOTS, :RibbonPolarization)  
#
#mt = init(BMCMSOTS, :RibbonPolarization_vsX; X=:braiding_time)
#
#pdata = map([merge(t0.get_plotparams(P),Dict("k"=>2.3)),
#		 delete!(t0.get_plotparams(P),"k")]) do P1 
#
#mt.plot(t0.get_plotparams(P))
#
#end 
#
#
#
##myPlots.plot(tasks)
