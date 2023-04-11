import myLibs: ComputeTasks, Utils
#import myPlots 

import LinearAlgebra,Statistics,Optim

import BMCMSOTS.WLO: select_mesh_point

import BMCMSOTS: MB

import PyPlot 

import BMCMSOTS:WLO, CalcWLOadapt, ChecksWLO ,BBH


import .CalcWLOadapt: fill_gaps_ks!, sum_kSteps_dist2pi!, rescaling_factor_dk, sum_kSteps ,init_gaps_ks, verify_dk_bounds 


P = (
		 braiding_time = 0.25, 
#		 braiding_time = rand(), 
		 s0_Hamilt = 0.1,
		 s_Hamilt = 1,
		 b_Hamilt = 1,
		 kPoint_start = -1, 
		 kMesh_type="Adaptive",
		 kMesh_model="line",
		 nr_kPoints = 21,
		 preserved_symmetries = "All",
#		 nr_perturb_strength = 3, 
#		 max_perturb_strength = 0.4, 
#		 nr_perturb_instances = 1, 
		 perturb_strength = 0.2,
		 )





nk = CalcWLOadapt.nr_kPoints(P)
k0 = CalcWLOadapt.kPoint_start(P)

Hdata = MB.get_pertHdata(MB.params_fromP(P), MB.H);


WF_line = WLO.init_storage1(WLO.init_eigH(rand(2),Hdata...)[1], nk)

WF_mesh = WLO.init_storage(WLO.init_eigH(rand(2),Hdata...)[1], nk)


dir1 = 2 # WLO in x direction, spectrum dep on ky


dir2 = 3-dir1  # perp direction 



kij = WLO.get_kij(nk,k0) 

WLO.store_on_mesh!!(WLO.eigH!, nk, kij, WF_mesh, Hdata...)

wlo1_mesh = WLO.wlo1_on_mesh_inplace(dir1, WF_mesh) 

@testset "WLO on mesh or line" begin 
	
	for k2=1:nk-1 

		get_K_ = WLO.get_kij_constrained(nk, k0, 2, kij(k2)[1], dir2)

		WLO.store_on_mesh1!!(WLO.eigH!, nk, get_K_, WF_line, Hdata...)

		wlo1_line = WLO.wlo1_on_mesh1_inplace(WF_line) 
		
		for k1=1:nk-1
	
			@test select_mesh_point(wlo1_mesh, WLO.orderinds(dir1,k1,k2))≈select_mesh_point(wlo1_line, k1) 
	
	end 
	end 
	
end 
		

k20 = rand(1:nk-1)

get_K = WLO.get_kij_constrained(nk, k0, 2, kij(k20)[1], dir2)

WLO.store_on_mesh1!!(WLO.eigH!, nk, get_K, WF_line, Hdata...)

WF_line_occ = WLO.psi_sorted_energy(WF_line; halfspace=true, occupied=true)
	
overlaps = WLO.init_overlaps_line(WF_line_occ)
overlaps2 = WLO.init_overlaps_line(WF_line_occ)

W1_line = WLO.init_wlo_mesh(WF_line_occ)
W1_line2 = WLO.init_wlo_mesh(WF_line_occ)

@show size.(overlaps) size(W1_line)




WLO.wlo1_on_mesh1_inplace!(W1_line, overlaps..., WF_line_occ) 

for k=1:nk-1 

	copy!(select_mesh_point(W1_line2,k),
				WLO.wlo1_one_inplace!(overlaps2..., WF_line_occ, k)
				)

end 

@testset "WLO line: one or all" begin 

	for k=1:nk-1 
	
		@test select_mesh_point(overlaps[1],k)≈ select_mesh_point(overlaps2[1],k)
		@test select_mesh_point(W1_line,k)≈ select_mesh_point(W1_line2,k)

		@test sort(WLO.get_periodic_eigvals(select_mesh_point(W1_line,k))) ≈ sort!(WLO.get_periodic_eigvals(select_mesh_point(W1_line,1)))


	end 

end 





results = CalcWLOadapt.Compute(P; observables=input_checks[:observables])


#import PyPlot  
#
#PyPlot.close()
#
#for y in eachcol(results["WannierBands1"][:,:,1])
#       PyPlot.plot(results["ks"],y)
#end


nk,kij,data_gap, data_2,dk_from_gap_1= results["data"];

#error()

#gap_at_k(k) = CalcWLOadapt.calc_gap!(data, k)
gap_at_k = CalcWLOadapt.calc_gap! 

gaps = [CalcWLOadapt.calc_gap!(data_gap, k) for k=kij(1:nk-1)]


#for k in kij(1:nk)


#end 

#gaps[1]=1

X=[0.02,0.5]
Y=[1e-8,pi/10]

X1 = LinRange(extrema(gaps)...,200)



PyPlot.close()

fig,ax=PyPlot.subplots(2,2,sharex=true,sharey=true)

#xs = LinRange(1e-10,1e-8,500);
#xs = LinRange(1e-10,0.3,500);
#
#
#
#PyPlot.plot(xs,[CalcWLOadapt.hardbound(x, Y...) for x=xs],label="hard")
#PyPlot.plot(xs,[CalcWLOadapt.softbound(x, Y...) for x=xs],label="soft")
#
#PyPlot.legend()
#




pdata = map(["square","sin", "line","expminv"],ax) do n,a

#	n=="line"||return 0 


	getf = getproperty(CalcWLOadapt, Symbol("f_"*n))

	getp = getproperty(CalcWLOadapt, Symbol(n*"_par_from_vals"))

	dk_from_gap = getf(getp(X,Y))

#	PyPlot.plot(X1, dk_from_gap.(X1), label=n)



	N=2nk 

println()


	gaps_new, ks_new = CalcWLOadapt.find_rescaling_factor_dk( N, k0, 
																													 (gap_at_k,data_gap), dk_from_gap , Y, )[1] 


	a.plot(ks_new, gaps_new, label=n, marker="o", markersize=2)

#	@show alpha 
#
##	PyPlot.plot(X1,f3.(X1),label=n)
##
##	PyPlot.scatter(gaps,f3.(gaps))
##



	a.legend()


end  
#
#for y in Y
#	PyPlot.plot(X,[y,y],c="gray")
#end 
#

#PyPlot.xlim(extrema(gaps))
#PyPlot.ylim(Y)















