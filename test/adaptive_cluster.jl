import myLibs: ComputeTasks, Utils
import myPlots 

import LinearAlgebra,Statistics,Optim

import BMCMSOTS.WLO: select_mesh_point

import BMCMSOTS: MB

import PyPlot 

import BMCMSOTS:WLO, CalcWLO 



#import .CalcWLO: fill_gaps_ks!, sum_kSteps_dist2pi!, sum_kSteps ,init_gaps_ks, verify_dk_bounds 



function p2D(kx_,ky_; kwargs...)


	PyPlot.close() 

	fig,ax = PyPlot.subplots(1)
	
	ax.set_xlim(-pi,pi)
	ax.set_ylim(-pi,pi)
	ax.set_aspect(1)
	fig.tight_layout()
	
	for kx in kx_
	
		ax.scatter(fill(kx,length(ky_)), ky_; kwargs...)
	
		#sleep(0.001) 
	
	end 
	
	
	ax.set_xlim(-pi,pi)
	ax.set_ylim(-pi,pi)
	ax.set_aspect(1)
	fig.tight_layout() 

end 

function p2D(kx,ky,nmax::Int; kwargs...)

	p2D(
			length(kx)<2nmax ? kx : kx[1:div(length(kx),nmax):end],
			length(ky)<2nmax ? ky : ky[1:div(length(ky),nmax):end];
			kwargs...
			)

end 


P = (braiding_time = 0.25, 
		 s0_Hamilt = 0.1, 
		 s_Hamilt = 1.0, #0.5, 
		 b_Hamilt = 1, 
		 nr_kPoints = 15, #1501
		 kMesh_model = "Line", 
		 kPoint_start = -1, 
		 preserved_symmetries = "All", 
		 perturb_strength = 0.0
		 )


plkw = (c="Olive",
				s=1,
				) 

kxy= CalcWLO.calc_kxy_adaptive(P) 

p2D(kxy...; plkw...)




































































































