using Distributed 
import LinearAlgebra 

@everywhere function H((kx,ky)::AbstractVector{<:Real})

	intra = [0.902 + 0.879im 0.391 + 0.623im 0.761 + 0.283im 0.277 + 0.26im; 0.376 + 0.355im 0.946 + 0.354im 0.95 + 0.834im 0.88 + 0.583im; 0.706 + 0.241im 0.782 + 0.785im 0.161 + 0.928im 0.023 + 0.559im; 0.495 + 0.933im 0.357 + 0.076im 0.728 + 0.825im 0.873 + 0.4im] 
	inter_x = [0.935 + 0.807im 0.821 + 0.786im 0.193 + 0.631im 0.983 + 0.415im; 0.624 + 0.878im 0.531 + 0.596im 0.262 + 0.108im 0.322 + 0.023im; 0.379 + 0.303im 0.094 + 0.479im 0.894 + 0.136im 0.206 + 0.992im; 0.071 + 0.624im 0.087 + 0.871im 0.938 + 0.819im 0.687 + 0.422im] 
	inter_y = [0.795 + 0.063im 0.014 + 0.86im 0.746 + 0.544im 0.349 + 0.185im; 0.912 + 0.427im 0.19 + 0.882im 0.032 + 0.609im 0.48 + 0.974im; 0.165 + 0.523im 0.201 + 0.424im 0.504 + 0.787im 0.538 + 0.583im; 0.16 + 0.363im 0.259 + 0.906im 0.056 + 0.317im 0.657 + 0.245im]

	LinearAlgebra.Hermitian(intra+cis(kx)*inter_x+cis(ky)*inter_y)

end 



pmap(H,eachcol(rand(2,10)))



