module WLO 
#############################################################################

import LinearAlgebra, Statistics 

import myLibs:Utils 

import ..Helpers: PeriodicFuns 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

# algo: https://www.johndcook.com/blog/standard_deviation/
# 1962 Welford  
# #
# used in: SnakeStates/src/DeviceWithLeads/src/TasksPlots.jl
# and in: Luxemburg/check_WLO/WLO/src/WLO.jl 
# 
# to be moved to myLibs!
#

function init_run_ms()::Vector{Float64}

	zeros(4)

end 

function init_run_m()::Vector{Float64}

	zeros(2)

end 

function init_run_ms!(storage::AbstractVector{Float64})::Nothing 

	storage .= 0 

	return 

end  

function run_m1!(storage::AbstractArray{<:Number}, 
								 k::Real, xk::Number,
								 inds...)::Nothing 

	storage[inds...] += (xk-storage[inds...])/k 

	return 

end 


function run_m!(storage::AbstractVector{Float64}, xk::Real)::Nothing

	storage[1] += 1

#	storage[2] += (xk-storage[2])/storage[1]

	run_m1!(storage, storage[1], xk, 2)

end 


function run_ms!(storage::AbstractVector{Float64}, xk::Real)::Nothing

	# storage[1]: k-1 
	
	storage[1] += 1

	# storage[1]: k

	#	storage[2] and storage[3]: M(k-1)

	storage[3] += (xk-storage[3])/storage[1]

	# storage[3]: M(k) 

	#	storage[4] S(k-1) 
	
	storage[4] += (xk-storage[2])*(xk-storage[3])

	# storage[4]: S(k) 
	
	storage[2] = storage[3] 

	# storage[2]: M(k)
	
	return 

end 

function run_var(storage::AbstractVector{Float64})::Float64 

	storage[1]>1 ? storage[4]/(storage[1]-1) : 0.0

end 


function run_mean(storage::AbstractVector{Float64})::Float64 

	storage[2] 

end 

function run_sigma(storage::AbstractVector{Float64})::Float64 
	
	sqrt(run_var(storage))

end  


function run_ms!(storage::AbstractVector{Float64},
								 X::AbstractVector{<:Real})::Nothing 

	for x in X 

		run_ms!(storage, x) 

	end  

end  


function run_ms(X::AbstractVector{<:Real})::Vector{Float64}

	storage = init_run_ms()

	run_ms!(storage, X)

	return storage

end 










#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function wcc_stat!(storage::AbstractVector{Float64},
									 S::AbstractVector{Float64},
									 X::AbstractVector{<:Real},
									 quantized_values::AbstractVector{<:Real},
									 )::Nothing 



	init_run_ms!(storage) 

	setindex!(X, PeriodicFuns.closest_periodic_a(X[1], quantized_values, 1), 1)

	run_m!(storage, X[1])
	
	for i = 2:lastindex(X)

		setindex!(X, PeriodicFuns.closest_periodic_a(X[i], X[i-1], 1), i)

		run_m!(storage, X[i])

	end 


	setindex!(S, 
						PeriodicFuns.closest_periodic_a(run_mean(storage), quantized_values, 1), 
						1)

	setindex!(S, 
						Statistics.stdm(X, PeriodicFuns.closest_periodic_b(S[1], quantized_values, 1)),
						2)


#	setindex!(S, run_mean(storage), 1)

#	setindex!(S, run_sigma(storage), 2)

	return 

end 


function wcc_stat!(X::AbstractArray{<:Real,N},
									 quantized_values::AbstractVector{<:Real};
									dim::Int=1 
								 )::Array{Float64,N} where N

	# take vectors along axis 'dim' and interate the other axes 

	storage = init_run_m()  # same for all sets of WCCs 


	S = zeros((size(X)[1:dim-1]..., 2, size(X)[dim+1:N]...))

	for I2 in CartesianIndices(axes(X)[dim+1:N])

		for I1 in CartesianIndices(axes(X)[1:dim-1])

			wcc_stat!(storage, 
								view(S, I1, :, I2), 
								view(X, I1, :, I2),
								quantized_values)

		end 

	end 

	return S

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




#function matmul!(A::AbstractMatrix{T},
#								 B::AbstractMatrix{Tb},
#								 C::AbstractMatrix{Tc},
#								 f1::Function=identity
#								 )::Nothing where {T<:Number,Tb<:T,Tc<:T}
#
#	A .= 0 
#
#	for k in axes(C,2), j in axes(B,2), i in axes(B,1)
#
#		A[i,k] += f1(B[i,j]) * C[j,k]
#
#	end 
#
#	return 
#
#end  
#

function overlap(wf1::AbstractMatrix{T1},
								 wf2::AbstractMatrix{T2}=wf1;
								 dim::Int=2
								 )::Matrix{promote_type(T1,T2)
													 } where {T1<:Number,T2<:Number} 

	A = Matrix{promote_type(T1,T2)}(undef, size(wf1,dim), size(wf2,dim))

	overlap!(A,wf1,wf2; dim=dim)

	return A 

end 



function overlap!(A::AbstractMatrix{<:Number},
									wf1::AbstractMatrix{<:Number},
								 wf2::AbstractMatrix{<:Number}=wf1;
								 dim::Int=2
								 )::Nothing 

	if dim==2 
		
		LinearAlgebra.mul!(A, wf1', wf2)

	elseif dim==1 
		
		LinearAlgebra.mul!(A, wf1, wf2')

		conj!(A)

	end 


	return 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function successive_overlaps(psi_W1::AbstractArray{T,2};
														 dim=nothing,dim_cases=nothing,
											 )::Array{T,3} where T<:Number  

	successive_overlaps_(psi_W1, dim, dim_cases)

end 


function successive_overlaps_(psi_W1::AbstractMatrix{T},
															dim::Nothing, dim_cases::Int
															)::Array{T,3} where T<:Number 

	# WFs have just one component --> 3rd (trivial) axis 
	
	@assert dim_cases in 1:2 

	successive_overlaps(reshape(psi_W1, size(psi_W1)..., 1); 
											dim_cases=dim_cases, dim=3)

end 



function successive_overlaps_(psi_W1::AbstractMatrix{T},
															dim::Int, dim_cases::Nothing
															)::Array{T,3} where T<:Number

	successive_overlaps_(psi_W1, dim, 3-dim)

end 


function successive_overlaps_(psi_W1::AbstractMatrix{T},
															dim::Int, dim_cases::Int
															)::Array{T,3} where T<:Number

	# there is just one WF with several components on axis 'dim'
	
	@assert dim in 1:2 
	@assert dim_cases==3-dim 

	successive_overlaps(reshape(psi_W1, size(psi_W1)..., 1);
												dim=dim,
												dim_cases=3-dim, 
												)

end 




function successive_overlaps(psi_W1::AbstractArray{T,3};
											 dim::Int,
											 dim_cases::Int,
											 repeat_first::Bool=true,
											 )::Array{T,3} where T<:Number 

	@assert dim!=dim_cases 

	overlaps_W1 = Array{T,3}(undef, 
													 size(psi_W1, dim), size(psi_W1, dim),
													 size(psi_W1, dim_cases) - 1 + repeat_first)
	


	for j in 1:size(psi_W1, dim_cases)-1 
	
		overlap!(selectdim(overlaps_W1, 3, j),
								 selectdim(psi_W1, dim_cases, j),
								 selectdim(psi_W1, dim_cases, j+1);
								 dim=dim-(dim_cases<dim))
	end 

	if repeat_first

		overlap!(selectdim(overlaps_W1,3,size(overlaps_W1,3)),
								 selectdim(psi_W1, dim_cases, size(psi_W1,dim_cases)),
								 selectdim(psi_W1, dim_cases, 1);
								 dim=dim-(dim_cases<dim))
	end 

	return overlaps_W1 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function occupied_subspace(k::AbstractVector, H::Function, data...;
													 kwargs...)#::Tuple{Matrix,Vector}

	occupied_subspace(H(k, data...); kwargs...)

end  


function init_occup(H::AbstractMatrix{T}
									 ) where T<:Number 
#									 )::Matrix{T} where T<:Number 

	n = LinearAlgebra.checksquare(H)

	return [Matrix{promote_type(T,Float64)}(undef, n, div(n,2)), 
					Vector{Float64}(undef, div(n,2))]

end 


function occupied_subspace(M::AbstractMatrix{T}, args...; kwargs...
														) where T<:Number 
#														)::Matrix{T} where T<:Number 

	u = init_occup(M)

	occupied_subspace!(u, M, args...; kwargs...)

	return u 

end   

function unoccupied_subspace(M::AbstractMatrix{T}, args...; kwargs...
														) where T<:Number 
#														)::Matrix{T} where T<:Number 

	u = init_occup(M)

	unoccupied_subspace!(u, M, args...; kwargs...)

	return u 

end   





function occupied_subspace!(u,
														#u::AbstractMatrix,
														k::AbstractVector, 
														H::Function, data...;
														kwargs...)::Nothing 

	occupied_subspace!(u, H(k, data...); kwargs...)

end  



function occupied_subspace!(u,#::AbstractMatrix,
														H::AbstractMatrix;
														kwargs...
														)::Nothing 

	eig = LinearAlgebra.eigen(H)

	occupied_subspace!(u, eig.vectors, eig.values; kwargs...)

end  


function unoccupied_subspace!(u,#::AbstractMatrix,
														H::AbstractMatrix;
														kwargs...
														)::Nothing 

	eig = LinearAlgebra.eigen(H)

	unoccupied_subspace!(u, eig.vectors, eig.values; kwargs...)

end  

function unoccupied_subspace!((u,v)::AbstractVector{<:AbstractArray},
#														u,::AbstractMatrix,
														psi::AbstractMatrix,
														E::AbstractVector{<:Real};
														kwargs...
														)::Nothing 

i = partialsortperm(E, div(length(E),2)+1:length(E); kwargs...)


	setindex!(u, selectdim(psi, 2, i), :, :)

	setindex!(v, view(E, i), :)


#	@show first(first(uv)) 
	
#	@show first(last(uv)) 


	return 

end 



function occupied_subspace!(uv::AbstractVector{<:AbstractArray},
#														u,::AbstractMatrix,
														psi::AbstractMatrix,
														E::AbstractVector{<:Real};
														kwargs...
														)::Nothing 

#println("here3")

#	@show first(first(uv))
#	@show first(last(uv)) 

	u,v = uv  

#	@show size(u) size(v)
#
#	@show E 


	i = partialsortperm(E, 1:div(length(E),2); kwargs...)

#	@show i kwargs


	setindex!(u, selectdim(psi, 2, i), :, :)

	setindex!(v, view(E, i), :)


#	@show first(first(uv)) 
	
#	@show first(last(uv)) 


	return 

end 

function occupied_subspace!(
														u::AbstractMatrix,
														psi::AbstractMatrix,
														E::AbstractVector{<:Real};
														kwargs...
														)::Nothing 

	occupied_subspace!([u,zeros(size(u,2))], psi, E; kwargs...)

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function unitary_part(G)

	F = copy(G) 

	unitary_part!(F) 

	return F 

end  

function unitary_part!(G)

#	return 

	svd = LinearAlgebra.svd!(G)  

	LinearAlgebra.mul!(G, svd.U, svd.Vt)

end 

function unitary_overlap(wf1::AbstractMatrix{T1},
								 wf2::AbstractMatrix{T2}
								 )::Matrix{promote_type(T1,T2)
													 } where {T1<:Number,T2<:Number} 

	A = overlap(wf1,wf2)

	unitary_part!(A)

	return A 

end 

function unitary_overlap!(A::AbstractMatrix{T},
									wf1::AbstractMatrix{T1},
								 wf2::AbstractMatrix{T2}
								 )::Nothing where {T<:Number,T1<:T,T2<:T}

	overlap!(A, wf1,wf2)

	unitary_part!(A) 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function check_start_int(n::Int, k0::Real)

	m0 = k0*(n-1)/pi 

	if !isapprox(m0,Int(round(m0)),atol=1e-10)

		@warn "Not all k-s will have mirror images in the mesh"

	end  

end 

function check_nr_kPoints(n::Int)

	@assert n>2 "Too few k points" 

end 


function get_kij(n::Int, k0::Real)::Function 

	check_nr_kPoints(n)

	check_start_int(n, k0)


	get_kij_((i,j)::NTuple{2,Int})::Vector{Float64} = get_kij_(i,j)

	function get_kij_(i::Int, j::Int)::Vector{Float64}
	
		@assert 1 <= i < n
		@assert 1 <= j < n
		
		k = [Float64(i-1), Float64(j-1)]
		
		k .*= -2pi/(n-1) 
	
		k .+= 2pi + k0
	
		return k
	
	end  

	return get_kij_

end 

function ind_minusk(i::Int, n::Int, k0::Real)::Int 

	PeriodicFuns.bring_periodic_to_interval(Int(round(k0*(n-1)/pi))+2-i, 1, n-1, n-1)

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function select_mesh_point(data::AbstractArray{T,M},
													 i::Int,
													 j::Int,
													 )::AbstractArray{T,M-2} where {T<:Number,M}
	
	selectdim(selectdim(data,M,j),M-1,i)

end 

function get_one_wrapper(get_one::Function,
												 data::AbstractArray{<:Number,M},
												 i::Int, 
												 j::Int,
												 data_...
												) where M 

	get_one(select_mesh_point(data, i, j), data_...)

end 

function get_one_wrapper(get_one::Function,
												 data::NTuple{N, <:AbstractArray} where N,
												 i::Int, 
												 j::Int,
												 data_...
												) 

	get_one(Tuple(select_mesh_point(A, i, j) for A in data), data_...)

end 




function get_one_wrapper(get_one::Function,
												 ij_to_arg::Function,
												 i::Int, 
												 j::Int,
												 data_...
												)

	get_one(ij_to_arg(i,j), data_...)

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function init_storage(item1::AbstractArray{T,N}, n::Int
											)::Array{T,N+2} where T<:Number where N

	repeat(item1, ones(Int, N)..., n-1, n-1)

end 

function init_storage(item1::Tuple, args...)::Vector

	[init_storage(it, args...) for it in item1]

end 


function store_on_mesh(get_one::Function, 
												arg2::AbstractArray{<:Number},
											 data...
												)


	store_on_mesh(get_one, nr_kPoints_from_mesh(arg2), arg2, data...)

end   


function store_on_mesh(get_one::Function, 
											 n::Int,
												arg2::Union{Function, 
																		AbstractArray{<:Number},
																		NTuple{N, <:AbstractArray} where N,
																		},
											 data...
												)


	storage = init_storage(get_one_wrapper(get_one, arg2, 1,1, data...), n)



	store_on_mesh!(get_one, n, arg2, storage, data...)

	return storage

end  



function store_on_mesh!(get_one::Function, 
												n::Int,
												arg2::Union{Function, AbstractArray{<:Number},
																		NTuple{N, <:AbstractArray} where N,
																		},
												data...
												)::Nothing 

	for j=1:n-1, i=1:n-1

		store_on_mesh_one!(get_one, arg2, i, j, data...)

	end 

	return 

end 






function store_on_mesh_one!(get_one::Function, 
												arg2::Union{Function, AbstractArray{<:Number},
																		NTuple{M, <:AbstractArray} where M,
																		},
														i::Int, j::Int,
												storage::AbstractArray{<:Number,N},
												data...
												)::Nothing where N

	setindex!(storage, 
						get_one_wrapper(get_one, arg2, i, j, data...),
						axes(storage)[1:N-2]..., i, j)

	return 

end  



function store_on_mesh_one!(get_one::Function, 
												arg2::Union{Function, AbstractArray{<:Number}},
														i::Int, j::Int,
														storage::AbstractArray{<:AbstractArray},
												data...
												)::Nothing 

	for (k,item) in enumerate(get_one_wrapper(get_one, arg2, i, j, data...))

		setindex!(storage[k], item, axes(storage[k])[1:end-2]..., i, j)

	end 

	return 

end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function get_item_ij((i,j)::Union{AbstractVector{Int},
																	Tuple{Int,Int}},
										 storage::AbstractArray{T,N}
										 )::AbstractArray{T,N-2} where T<:Number where N

	get_item_ij(i,j, storage)

end 

function get_item_ij(i::Int, j::Int, 
										 storage::AbstractArray{T,N}
										 )::AbstractArray{T,N-2} where T<:Number where N

	PeriodicFuns.cycseldim(PeriodicFuns.cycseldim(storage, N, j), N-1, i)

end 

function get_item_ij(k::Int, dir::Int,
										 start_ij::Union{AbstractVector{Int},NTuple{2,Int}},
										 storage::AbstractArray{T,N}
										 )::AbstractArray{T,N-2} where T<:Number where N

	get_item_ij([d==dir ? s+k-1 : s for (d,s)=enumerate(start_ij)], storage)

end 

function get_item_line(dir::Int, k::Int, 
										 storage::AbstractArray{T,N}
										 )::AbstractArray{T,N-1} where T<:Number where N

	selectdim(storage, N-2+dir, k)

end  # k[dir] 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#
#end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function wlo(get_wf::Function,
							 get_wf!::Function,
							 k_in::AbstractVector{<:Real},
							dir::Int,
						 data...)
	
	k_fin = copy(k_in)

	k_fin[dir] += 2pi 
	
	return wlo(get_wf, get_wf!, k_in, k_fin, data...)

end 


function wlo(get_wf::Function,
							 get_wf!::Function,
							 k::Real,
							dir::Int,
						 args...)

	wlo(get_wf, get_wf!, construct_k(k, dir), dir, args...)

end 



function wlo1(args...)

	wlo(first∘occupied_subspace, occupied_subspace!, args...)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function wlo_from_mesh(
											 IJ0::Union{AbstractVector{Int},Tuple{Int,Int}},
											 dir::Int,
											 wfs::AbstractArray{ComplexF64,4},
											 )::Matrix{ComplexF64} 

	@assert 1<=dir<=2 "Direction must be 1=='x' or 2=='y'"

	return mapreduce(*, 1:nr_kPoints_from_mesh(wfs)-1) do k 

		unitary_overlap(get_item_ij(k, dir, IJ0, wfs),
										get_item_ij(k+1, dir, IJ0, wfs))


	end 
	
end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function wlo1_on_mesh(dir1::Int, n::Int, k0::Real, Hdata...)

	Bloch_WFs = store_on_mesh(first∘occupied_subspace, n, 
														get_kij(n,k0), Hdata...)

	return store_on_mesh(wlo_from_mesh, tuple, dir1, Bloch_WFs), Bloch_WFs 

end 


function wlo1_on_mesh(dir1::AbstractVector{Int}, 
											n::Int, k0::Real,
											Hdata...)

	Bloch_WFs = store_on_mesh(first∘occupied_subspace, n, 
														get_kij(n,k0), Hdata...)


	w1 = store_on_mesh(wlo_from_mesh, n, tuple, dir1[1], Bloch_WFs)


	return (
					(w1, (wlo1_on_mesh(d1, n, Bloch_WFs) for d1 in dir1[2:end])...), 
					Bloch_WFs
					)

end 



function wlo1_on_mesh(dir1::Int, n::Int, Bloch_WFs::AbstractArray,
											)

	store_on_mesh(wlo_from_mesh, n, tuple, dir1, Bloch_WFs)


end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function Wannier_subspace_on_mesh_one(ij::NTuple{2,Int},
																			W::AbstractArray,
													s::Union{Function,Real}, 
													)

	Wannier_subspace(get_item_ij(ij, W), s)

end 

	
function nr_kPoints_from_mesh(W::AbstractArray)::Int 
	
	only(unique(size(W)[end-1:end]))+1

end 


function Wannier_subspace_on_mesh(W::AbstractArray,
													s::Union{Function,Real}, 
													)

	store_on_mesh(Wannier_subspace_on_mesh_one, 
								nr_kPoints_from_mesh(W), tuple, W, s)

end 


function Wannier_subspace_on_mesh(s::Union{Function,Real}, Hdata...)

	W,U = wlo1_on_mesh(Hdata...)

	return Wannier_subspace_on_mesh(W, s),W,U

end 


#function Wannier_band_basis_mesh(Bloch_WFs::AbstractArray,
#																 WLO_WFs::AbstractArray)
#
#	store_on_mesh(Wannier_band_basis0, tuple, Bloch_WFs, WLO_WFs)
#
#end  

function Wannier_band_basis_mesh(
																 W::AbstractArray,
																 Bloch_WFs::AbstractArray,
																 s::Union{Function,Real})
																 
	eig = Wannier_subspace_on_mesh(W, s)

	return Wannier_band_basis_mesh(eig[1], Bloch_WFs)[1], eig 

end  

function Wannier_band_basis_mesh(
																 eigvecs::AbstractArray,
																 Bloch_WFs::AbstractArray
																 )
	(															 
	store_on_mesh(Wannier_band_basis0, 
								nr_kPoints_from_mesh(Bloch_WFs), 
								tuple, Bloch_WFs, eigvecs),
	eigvecs
	)

end  



function Wannier_band_basis_mesh(
																 s::Union{Function,Real},
																 data...)

	Wannier_band_basis_mesh(wlo1_on_mesh(data...)..., s)

end  

	
function wlo2_on_mesh(					 
											dir2::Int,
											data...)

	nested_WF_dir1 = Wannier_band_basis_mesh(data...)[1]

	return store_on_mesh(wlo_from_mesh, 
											 nr_kPoints_from_mesh(nested_WF_dir1),
											 tuple, dir2, nested_WF_dir1)
	
end 


function wlo2_on_mesh(dir2::Int,
											sector::Union{Function,Real},
											H::Function,
											Hdata...)

	wlo2_on_mesh(dir2, sector, [2,1][dir2], H, Hdata...)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function wlo(get_wf::Function,
							 get_wf!::Function,
							k_in::AbstractVector, 
							k_fin::AbstractVector, data...)


#	pr(args...) = [data[1]==1 && println(string(args...))]



	k, kdata = prep_feed_ks(k_in, k_fin)  


	u_old = get_wf(k, data...)

#	pr("u_old and u1 at k=",k)


	update_k!(k, kdata) 

	u_new = get_wf(k, data...)
	
#	pr("u_new at k=",k)



	ov = unitary_overlap(u_old, u_new)
#	if data[1]==1 
#
#		ov[2,1]=0
#		ov[1,2]=0
#	end 

	W = copy(ov) 

	#pr("size(W)=",size(W))

	u1 = copy(u_old)

	s = sum(u1)
#storage: 	u_old u_new u1 ov W 

	while update_k!(k, kdata)


		u_old .= u_new  

		get_wf!(u_new, k, data...)

		unitary_overlap!(ov, u_old, u_new) 

#	if data[1]==1 
#
#		ov[2,1]=0
#		ov[1,2]=0
#	end 

		W *= ov 

	end 
	
#	pr("last (unused) k: ",k)
	
#	@assert s≈ sum(u1) 


	unitary_overlap!(ov, u_new, u1)
#	if data[1]==1 
#
#		ov[2,1]=0
#		ov[1,2]=0
#	end 

	W *= ov 

#	pr("\n")

	return W 
	
end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Wannier_Hamilt(W::AbstractMatrix{<:Number})::Matrix{ComplexF64}

	H = log(Array{ComplexF64}(W))

	H *= -im 

	return H 

end 






function get_periodic_eigvals(W::AbstractMatrix)::Vector{Float64}

	get_periodic_eigvals(LinearAlgebra.eigvals(W))

end  

function get_periodic_eigvals(eigvals::AbstractVector)

	get_periodic_eigvals.(eigvals)

end 

function get_periodic_eigvals(eigval::Number)::Float64

	angle(eigval)/2pi

end 








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function Wannier_min_gap(W::AbstractMatrix)::Float64
	
	nu = get_periodic_eigvals(W)

	D = [Utils.dist_periodic(nu[i],nu[j],1) for i=1:length(nu) for j=i+1:length(nu)] 

	return minimum(D) 

end 

function Wannier_min_gap(args...)::Float64 

	Wannier_min_gap(wlo1(args...))

end  


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function choose_Wannier_subspace(s::Real)

	Dict(:rev=>s>0)

end 


function choose_Wannier_subspace(s::Function)

	choose_Wannier_subspace(s(1))

end 


function Wannier_subspace(W::AbstractMatrix,
													s::Union{Function,Real}, 
													)#::Matrix 

	eig = LinearAlgebra.eigen(W)


	ov = overlap(eig.vectors)

#	abs(ov[2,1]) + abs(ov[1,2]) <1e-12 || @warn "WLO eigvec not orthogonal"

#	round.(real(),digits=3)

#	@show size(eig.values) size(eig.vectors)

	@assert LinearAlgebra.norm(W * eig.vectors - eig.vectors * LinearAlgebra.Diagonal(eig.values)) < 1e-12 

	#@assert W * eig.vectors ≈ eig.vectors * LinearAlgebra.Diagonal(eig.values)
	


#	println("here1")

	occ, vocc = occupied_subspace(eig.vectors, 
														 get_periodic_eigvals(eig.values);
														 choose_Wannier_subspace(s)...)


#	@show size(occ) size(vocc) 
	
#	@show eig.values
#	@show get_periodic_eigvals(eig.values) vocc

#	n = sum(abs, eig.values)/2 

#	@show n 


	@assert any([LinearAlgebra.norm(W * occ -occ * LinearAlgebra.Diagonal(n*exp.(im*2*pi*vocc)))<1e-10 for n in abs.(eig.values)])
	
#	LinearAlgebra.norm(W * occ -occ * LinearAlgebra.Diagonal(n*exp.(im*2*pi*vocc)))<1e-10 || @warn ""

#	 W * occ ≈ occ * LinearAlgebra.Diagonal(n*exp.(im*2*pi*vocc))
#	@assert W * occ ≈ occ * LinearAlgebra.Diagonal(exp.(im*2*pi*vocc))


#	v = only(vocc)

#	@show v 
#
#	@show LinearAlgebra.norm(W*occ - exp(im*2*pi*v)*occ)
#	@show LinearAlgebra.norm(W*occ - exp(im*2*pi*(v+1))*occ)
#	@show LinearAlgebra.norm(W*occ - exp(im*2*pi*(v-1))*occ)


	if choose_Wannier_subspace(s)[:rev] 

		#only(vocc)>=0 || @warn ""

	else 

		#only(vocc)<=0 || @warn ""

	end 


	return occ, vocc

end 

















#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Wannier_band_basis0(ij::Tuple{Int,Int},
														 Bloch_WFs_::AbstractArray{ComplexF64,4},
														 W1wf_::AbstractArray{ComplexF64,4},
														 )::Matrix{ComplexF64}

	get_item_ij(ij, Bloch_WFs_) * get_item_ij(ij, W1wf_)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



# w_x in (6.6)

function Wannier_band_basis0(k::AbstractVector, 
														dir1::Int, 
														subspace::Union{Real,Function},
														H::Function,
														Hdata...)
	
	U,E = occupied_subspace(k, H, Hdata...)

#####
#	@show size(U)
#
#	@show abs.(U'*U)
#
#	GramSchmidt!(U)
#
#	@show abs.(U'*U) 
#
#	error()
#######

	w1 = wlo1(k, dir1, H, Hdata...) 

#######
#
#	W = LinearAlgebra.eigen(w1).vectors 
#
#	@show size(W)
#
#	@show round.(abs.(W'*W),digits=4)
#
##	GramSchmidt!(W,2)
#
##	@show abs.(W'*W) 

#	error()  
########


	W,nu = Wannier_subspace(w1, subspace) 

#####
#
#	@show size(W)
#
#	@show round.(abs.(W'*W),digits=4) 
#
#	error()
#
####

	out = U*W 

	return U,E,w1,W,nu,out

end 

function IPR(out::AbstractVector)::Float64

	(1/sum(abs2∘abs2, out) - 1)/(length(out) - 1)

end 


function Wannier_band_basis(k::AbstractVector{<:Real}, 
														dir1::Int, 
														subspace::Union{Real,Function},
														storage::AbstractVector,
														H::Function,
														Hdata...)::Matrix{ComplexF64}


	psi_H, E, W1, psi_W1, nu, out = Wannier_band_basis0(k, dir1, subspace, H, Hdata...)

####
	@assert nu ≈ get_periodic_eigvals(LinearAlgebra.diag(psi_W1' * W1 * psi_W1))

		eig = LinearAlgebra.eigen(W1)

	@assert only(nu)==get_periodic_eigvals(eig.values[argmax(abs.(overlap(eig.vectors,psi_W1))[:])])

@assert E ≈ LinearAlgebra.diag(psi_H' * H(k, Hdata...) * psi_H)

#	error()

#	println(abs(only(diff(E))))
###

#@show size(LinearAlgebra.eigen(W1).vectors)

#	ov = abs.(overlap(LinearAlgebra.eigen(W1).vectors))
	
#	@show round.(ov,digits=5)

#@show size(psi_W1)
#	@show k 


	dir2 = 3-dir1  


	push!(storage, [only(nu),
											
									copy(k),

									dir2, 

									only(eachcol(out)),

									E,

									only(eachcol(psi_W1)), 

									psi_H 
									
									]

				)


### 



return copy(out)

end 



function Wannier_band_basis!(w::AbstractMatrix,
														 k::AbstractVector, 
														dir1::Int, 
														subspace::Union{Real,Function},
														Hdata...)::Nothing 


	w .= Wannier_band_basis(k, dir1, subspace, Hdata...) 

#	U = occupied_subspace(k, Hdata...)  

#	W = Wannier_subspace(wlo1(k, dir1, Hdata...), subspace)

#	LinearAlgebra.mul!(w, U, W) 

	return 
end 




function wlo2(k::Union{AbstractVector,Real}, dir2::Int, 
							subspace::Union{Real,Function},
						 args...)

	wlo(Wannier_band_basis, 
			Wannier_band_basis!,
			k, dir2, [2,1][dir2], subspace, args...)


end 

function wlo2_pm_evals(k::Union{AbstractVector,Real}, 
											 dir2::Int, 
											 args...)::Matrix{Float64}

	hcat(get_periodic_eigvals(wlo2(k, dir2, +, args...)),
			 get_periodic_eigvals(wlo2(k, dir2, -, args...))
			 )

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function wlo1_evals_pm(args...)

	sort!(get_periodic_eigvals(wlo1(args...)), rev=true)

end 



#===========================================================================#
#
# Some plotting functions 
#
#---------------------------------------------------------------------------#


function scatter_yperiodic(ax, x, y; kwargs...)

	for (yr,kw) in prep_yperiodic(y; kwargs...)

		ax.scatter(x,yr; kw...)

	end 

end 

function plot_yperiodic(ax, x, y; kwargs...)

	for (yr,kw) in prep_yperiodic(y; kwargs...)

		ax.plot(x,yr; kw...)

	end 

end 

function plot_wcc(ax, xy::Int, W_, K_; s::Real=10, loc=nothing,kwargs...)

	W,k = parse_input_wcc(W_,K_)


	for pm=1:2
		
		scatter_yperiodic(ax, k[pm,:], W[pm,:]; s=[s,sqrt(s)][pm], c=["k","r"][pm], zorder=5, 
											label=string("\$","+-"[pm],"\$"),
											kwargs...)
		
	end 


	d1 = "xy"[xy]

	d2 = "xy"[3-xy] 


	ax.set_xlabel("\$k_$d2/\\pi\$")

	ax.set_ylabel("\$\\nu_$d1\$",rotation=0)


	@show size(W)

	error("dist_modulo,diff -> Utils.dist_periodic") 

#	x->Utils.dist_periodic(0,w,1)

	gap = max(1e-12, minimum(PeriodicFuns.dist_modulo, diff(W,dims=1)))


	gap = gap<1e-10 ? "0.0" : string(round(gap, digits = Int(ceil(-log10(gap)))+2))[1:min(end,6)]

	gap = rpad(gap, 6, "0")

	ax.set_title("\$W_$d1\$  (gap\$=$gap\$)")
	
	ax.set_xlim(k0/pi .+ [0,2])

	ax.set_ylim(-0.7, 0.7)

	if all(!isempty, ax.get_legend_handles_labels())

		ax.legend(loc=loc)

	end 

	return 

end  


function parse_input_wcc(W_,K_)
	
	parse_input_wcc(W_), parse_input_wcc(K_)/pi  

end 



function plot_wcc_norm(ax, xy::Int, W_,K_, k0; s::Real=10, loc=nothing,kwargs...)

	W,k = parse_input_wcc(W_,K_)


	for pm=1:2

		ax.scatter(k[pm,:], W[pm,:]; s=s, c=["k","r"][pm], zorder=5, 
							 label=string("\$|\\lambda_","+-"[pm],"|\$"),
											kwargs...)
		
	end 


	ax.set_xlim(k0/pi .+ [0,2])

	d1 = "xy"[xy]

	d2 = "xy"[3-xy] 

	ax.set_xlabel("\$k_$d2/\\pi\$")

	ax.set_ylabel("\$\\nu_$d1\$",rotation=0)

	ax.set_title("\$W_$d1\$: eigenvalue norm")
	
	ax.set_ylim(-.05,1.05)
	
	if all(!isempty, ax.get_legend_handles_labels())

		ax.legend(loc=loc)

	end  

	return 

end  


function parse_input_wcc(A::AbstractMatrix)

	copy(A)

end 

function parse_input_wcc(A::AbstractVector)

	vcat(transpose(A),transpose(A))

end  




function prep_yperiodic(y;
												yreps::AbstractVector{<:Real}=[-1,0,1],
													 alpha_max::Real=1, alpha_min::Real= 0.4,  
													 alpha=nothing, label=nothing, 
													 kwargs...)
	
	alpha_min = min(alpha_min,alpha_max)

	r_min,r_max = extrema(abs,yreps)
	
	a = (alpha_min-alpha_max)/(r_max-r_min)

	b = alpha_max - a*r_min 

	return [(y .+ rep,
					 merge(Dict(kwargs),
								 Dict(:label=>(rep==0 ? label : nothing),
											:alpha=>a*abs(rep)+b
											)
								 )
					 )  for rep in yreps]

end 

function plot_wcc2(ax, xy::Int, W_, K_, k0; s::Real=10,
									 loc=nothing,
									kwargs...)

	W,k= parse_input_wcc(W_,K_)


	klim = k0/pi .+ [0,2]


	for pm=1:2 
	
		scatter_yperiodic(ax, k[pm,:], W[pm,:]; s=[s,sqrt(s)][pm], c=["k","r"][pm], 
											label=string("\$","+-"[pm],"\$"),
											kwargs...)

	end 

	plot_yperiodic(ax, klim, fill(0.5,2);
								 color="k",alpha_max=0.2,zorder=get(kwargs,:zorder,0)-10)




	d1 = "xy"[xy]

	d2 = "xy"[3-xy] 

	ax.set_xlabel("\$k_$d1/\\pi\$")

	ax.set_ylabel("\$\\phi_{$d1$d2}\$",rotation=0) 
	
	ax.set_title("\$W_{$d1$d2}\$",rotation=0) 

	ax.set_xlim(klim)
	
	ax.set_ylim(-0.7, 0.7)


	if all(!isempty, ax.get_legend_handles_labels())

		ax.legend(loc=loc)

	end 


end 




















































































































































































































































































#############################################################################
end # module WLO 
