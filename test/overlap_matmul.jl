import BMCMSOTS: WLO 
import LinearAlgebra,PyPlot
import BMCMSOTS.TasksPlots: line, fit_line_conv ,line_intersection


@testset "multiply pairs recursively" begin 

	for N= 1:100#(10^1) .+ (1:10)#[10^6,1+10^6,2+10^6]#1:873:100000
	
		matrices = [rand(ComplexF64,3,3) for i=1:N]
		
		storage = rand(ComplexF64,3,3,2,N)
		
		for i=1:N
		
			copy!(view(storage,:,:,1,i),matrices[i])
		
		end 
		
		result1 = reduce(*,matrices)
		
		@test  result1 ≈ WLO.matmulpairs!(selectdim(storage,3,2),
		 																		selectdim(storage,3,1))
		
		
		src = rand(ComplexF64,3,3,N)
		dst = rand(ComplexF64,3,3,N)
		
		for i=1:N
		
			copy!(selectdim(src,3,i),matrices[i])
		
		end 
		
		A = copy(src)
		@test result1 ≈ WLO.matmulpairs!(dst,src)
		@test A≈src 
		
		@test result1≈WLO.matmulpairs!(dst, matrices...)


		src = rand(ComplexF64,3,3,N+1)

		for i=1:N
		
			copy!(selectdim(src,3,i),matrices[i])
		
		end 

		@test WLO.matmulpairs!(src,N)≈result1


	end 
	
end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



hsize = 32 
intra=rand(ComplexF64,hsize,hsize) 
inter_x=rand(ComplexF64,hsize,hsize) 
inter_y=rand(ComplexF64,hsize,hsize) 

function H((kx,ky)::AbstractVector{<:Real})

	LinearAlgebra.Hermitian(intra+cis(kx)*inter_x+cis(ky)*inter_y)

end 




#nmeshes = 10:20:100
nmeshes = Vector{Int}(Utils.uniqlogsp(4, 1000, 50, 0; Trunc=true))[1:20]

t1 = zeros(length(nmeshes))
t2 = zeros(length(nmeshes))
PyPlot.close()

@testset "overlap mult. methods" begin 

	for (i_n,nmesh) in enumerate(nmeshes)

		println()
		@show nmesh 

		psiH = WLO.psiH_on_mesh(nmesh,0, H) 
		
		dest1 = WLO.init_wlo_mesh(psiH,nmesh)
		dest2 = WLO.init_wlo_mesh(psiH,nmesh)
		
		overlaps = WLO.init_overlaps_line(psiH, nmesh) 
 
		k2 = rand(1:nmesh-1)

		dir1= rand(1:2)


#		WLO.unitary_overlaps_on_mesh!(overlaps[1], dir1, WLO.orderinds(dir1,1,k2), psiH)

		WLO.wlo1_on_mesh1_inplace!(dest1, overlaps..., psiH, dir1, k2)
		WLO.wlo1_on_mesh1_inplace_2!(dest2, overlaps..., psiH, dir1, k2)

		overlaps[1].=0
		overlaps[2].=0

		t1[i_n] = @elapsed WLO.wlo1_on_mesh1_inplace!(dest1, overlaps..., psiH, dir1, k2)
		overlaps[1].=0
		overlaps[2].=0
		t2[i_n] = @elapsed WLO.wlo1_on_mesh1_inplace_2!(dest2, overlaps..., psiH, dir1, k2)


		PyPlot.scatter(nmesh,t1[i_n],c="b")
		PyPlot.scatter(nmesh,t2[i_n],c="r")
		
		i_n==1 && PyPlot.legend() 

		@test dest1 ≈ dest2 rtol=1e-10 atol=1e-10

	end 	
	
	
	PyPlot.loglog(nmeshes,t1,c="b",label=1)
	PyPlot.loglog(nmeshes,t2,c="r",label=2)


	pars = map(enumerate((t1,t2))) do (i,t) 

		stop = findlast(>(1e-8),t)

		par,lab=fit_line_conv(log10.(nmeshes),log10.(t),1,stop)


		PyPlot.loglog(nmeshes, exp10.(line(log10.(nmeshes),par)), label="$i: $lab",ls="--")

		return par 

	end 

	n0, t0 = exp10.(line_intersection(pars...))


	PyPlot.scatter(n0,t0,zorder=10,c="pink",s=100,label="Intersection: n=$(Int(round(n0)))")


	PyPlot.legend() 

	PyPlot.xlabel("\$N\$")
	PyPlot.ylabel("\$t\$",rotation=0) 

	PyPlot.title("size(H)=$(size(H(rand(2))))")

	PyPlot.tight_layout()

	println() 

end  





