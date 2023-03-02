import BMCMSOTS: WLO 


@testset "multiply pairs recursively" begin 

	for N= 1:100#(10^1) .+ (1:10)#[10^6,1+10^6,2+10^6]#1:873:100000


matrices = [rand(3,3) for i=1:N]


storage = rand(3,3,2,N)


for i=1:N

	copy!(view(storage,:,:,1,i),matrices[i])

end 

result1 = reduce(*,matrices)

@test  result1 ≈ WLO.matmulpairs!(selectdim(storage,3,2),
 																		selectdim(storage,3,1))


src = rand(3,3,N)
dst = rand(3,3,N)

for i=1:N

	copy!(selectdim(src,3,i),matrices[i])

end 

A = copy(src)
@test result1 ≈ WLO.matmulpairs!(dst,src)
@test A≈src 

@test result1≈WLO.matmulpairs!(dst, matrices...)

end 

end 







