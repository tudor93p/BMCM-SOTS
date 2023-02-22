using BenchmarkTools 

import BMCMSOTS.Helpers: PeriodicFuns 

import myLibs:Utils 


L = 3


@testset "closest_periodic_shifted_a" begin 

	
	for i=-L:L 

		break 

		@show i
		
		for j=-L:L, k=0:2L+3
	
		for a in rand(5).-0.5, c in rand(5)

			T = max(k+c,1e-3) 

#			for b_ in rand(5)
			for b_ in [rand(5);collect(eachcol(rand(5,5)))]

				b = b_ .- 0.5  .+ j
	
			#		worst case scenario: T=1e-3, |a-b|=2L+1 => nmax = 2L+1/(1e-3)

			a1 = Utils.closest_periodic_shifted_a(i+a,b,T)
			
			 a2 = PeriodicFuns.closest_periodic_a2(i+a,b,T, 1000*(2L+5))

			@test isapprox(a1,a2,rtol=1e-5)|isapprox(a1,a2,atol=1e-5)

		end 
		end 
	
	end 
end 
	

end 

L = 7

@testset "closest_periodic_b" begin 

	
	for i=-L:L 

		@show i
		
		for j=-L:L, k=0:2L+3
	
		for a in rand(5).-0.5, c in rand(5)

			T = max(k+c,1e-3) 

			for b_ in collect(eachcol(rand(5,5)))

				b = j.+ b_ .- 0.5 
	
			#		worst case scenario: T=1e-3, |a-b|=2L+1 => nmax = 2L+1/(1e-3)

			b1 = Utils.closest_periodic_b(i+a,b,T)
			
			b2 = PeriodicFuns.closest_periodic_b2(i+a,b,T)

			@test isapprox(b1,b2,rtol=1e-5)|isapprox(b1,b2,atol=1e-5)

		end 

		end 
	
	end 
	
end 

end 

