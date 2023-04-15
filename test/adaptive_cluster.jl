import myLibs: ComputeTasks, Utils
import myPlots 

import LinearAlgebra,Statistics,Optim

import BMCMSOTS.WLO: select_mesh_point

import BMCMSOTS: MB

import PyPlot 

import BMCMSOTS:WLO, CalcWLO 



#import .CalcWLO: fill_gaps_ks!, sum_kSteps_dist2pi!, sum_kSteps ,init_gaps_ks, verify_dk_bounds 



function p2D((fig,ax), kx_,ky_; kwargs...)

	lim = [-1,1]

	
	ax.set_aspect(1)
	fig.tight_layout()
	
	for kx in kx_
	
		ax.scatter(fill(kx/pi,length(ky_)), ky_/pi; kwargs...)
	
		#sleep(0.001) 
	
	end 
	
	
	ax.set_xlim(lim)
	ax.set_ylim(lim)

	ax.set_xlabel("\$k_x\\;[\\pi]\$",labelpad=3)
	ax.set_ylabel("\$k_y\\;[\\pi]\$",labelpad=-2.6)#rotation=0)

	ax.set_aspect(1)
return fig,ax
end 

function p2D(arg,kx,ky,nmax::Int; kwargs...)

	p2D(arg,
			length(kx)<2nmax ? kx : kx[1:div(length(kx),nmax):end],
			length(ky)<2nmax ? ky : ky[1:div(length(ky),nmax):end];
			kwargs...
			)

end 


PyPlot.ioff()
#PyPlot.ion() 
PyPlot.close() 

fig,Ax = PyPlot.subplots(1,2,figsize=(8.6/2.54,0.53*8.6/2.54),
												 #sharex=true,sharey=true,
											 )

P = (braiding_time = 0.25, 
		 s0_Hamilt = 0.1, 
		 s_Hamilt = 1.0, #0.5, 
		 b_Hamilt = 1, 
		 nr_kPoints =81,# 15, #1501
		 kMesh_model = "Line", 
		 kPoint_start = -1, 
		 preserved_symmetries = "All", 
		 perturb_strength = 0.0
		 )


plkw = (c="Olive",
				s=1,
				) 

kxy= CalcWLO.calc_kxy_adaptive(P) 

for ax in Ax 
p2D((fig,ax), kxy...; plkw...)
end 

a2 = Ax[2]
a2.set_ylabel(nothing)
a2.set_yticklabels(["" for t=a2.get_yticks()])

for (g,s) in ((a2.get_ylim,a2.set_ylim),
							(a2.get_xlim,a2.set_xlim)
							)
	m,M = g() 
	s((m+M)/2,M)
end 

Ax[1].tick_params("y",pad=2) 

for ax in Ax 
for a in [ax.xaxis.label, ax.yaxis.label]
	a.set_fontsize(9)
 end 
 for A in (ax.get_xticklabels(),ax.get_yticklabels()), a in A 
	 a.set_fontsize(7.5)
	end 
end 

fig.subplots_adjust(
top=0.985,
bottom=0.18,
left=0.128,
right=0.976,
hspace=0.2,
wspace=0.105
)

fig.savefig("adaptive_mesh.pdf",format="pdf")


































































































