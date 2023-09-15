



import BMCMSOTS: ChiralPwave, CalcWLO



ChiralPwave.H(rand(2),rand(2),ChiralPwave.dvector_chiral)
ChiralPwave.H(rand(2),rand(2),ChiralPwave.dvector_helical)
ChiralPwave.H(rand(2),rand(2),ChiralPwave.psi_extswave)


obs = ["WannierBands1"] 



P = Dict(
				 :nr_kPoints =>30,

				 :Delta0=>0.3, 

#				 :pairSymm=>"ext.swave",
				 :pairSymm=>"chiral",

				 :pairSymmConfig=>1,


				 )




results = CalcWLO.Compute(P; observables=obs)

@testset "p-wave twofold degeneracy" begin 

	@test ≈(eachslice(results["WannierBands1"],dims=3)...) atol=1e-7

end 

