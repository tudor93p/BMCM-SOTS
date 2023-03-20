import myPlots 
import BMCMSOTS: RibbonWLO , MB , WLO

import myLibs:TBmodel,Lattices 

P = (width=10,
			braiding_time = 0.25, 
			s0_Hamilt = 0.1, s_Hamilt = 1, b_Hamilt = 1, 
			kPoint_start = -1, nr_kPoints = 10, 
			preserved_symmetries = "All",#Ct",#"None",#"Ct", 
			nr_perturb_strength = 3, max_perturb_strength = 0.4, 
			nr_perturb_instances = 1, perturb_strength = 0.2,) 

observables=RibbonWLO.calc_observables 

data = RibbonWLO.Compute(P; observables=observables)

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

t = init(BMCMSOTS, :RibbonPolarization) 

pdata = t.plot(t.get_plotparams(P)) 

t2 = init(BMCMSOTS, :RibbonWannierBands1); 

pdata2 = t2.plot(t2.get_plotparams(P)) 

t3 = init(BMCMSOTS, :RibbonWannierDensity);

pdata3 = t3.plot(t3.get_plotparams(P)) 


myPlots.plot(t2, t3, t)
