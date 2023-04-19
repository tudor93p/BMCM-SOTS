module WLO 
#############################################################################

import LinearAlgebra, Statistics 

import SharedArrays: SharedArray  

import Distributed: @spawnat, myid, workers,nworkers,procs
import DistributedArrays: 	DArray, SubOrDArray, localindices,localpart, dzeros 




import myLibs:Utils,SignalProcessing 
import myLibs.Parameters:UODict


import ..Helpers: PeriodicFuns, Symmetries 

const SubOrArray{T,N} = Union{Array{T,N}, SubArray{T,N,<:Array}}
const SubOrSArray{T,N} = Union{SharedArray{T,N}, SubArray{T,N,<:SharedArray}}

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

usedkeys::Vector{Symbol} = [:nr_kPoints, :kPoint_start]

function nr_kPoints(P::UODict)::Int 

	P[:nr_kPoints]

end  

function kPoint_start(P::UODict)::Float64

	pi*P[:kPoint_start]

end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function same_procs(ads::Tuple{Vararg{<:AbstractDict}})::Bool

	s1 = Set(keys(ads[1]))

	for i=2:length(ads) 

		s1==Set(keys(ads[i])) || return false 

	end  

	return true 

end 


function same_mesh_distrib(a::AbstractDict,
								b::AbstractDict
								)::Bool  

	same_procs((a,b)) || return false 

	for (k,va) in pairs(a) 

		va==b[k] || return false 

	end 

	return true 

end 



function same_mesh_distrib(A::SubOrDArray,
													 B::SubOrDArray,
													 meshdim::Int=2 
													 )::Bool 

	for (a,b) in zip(A.indices,B.indices), d=1:meshdim 

		select_mesh_dir(a,d,meshdim)==select_mesh_dir(b,d,meshdim) || return false 

	end 

	return true 

end  
	
function same_mesh_distrib(arrays::Union{<:AbstractVector{<:Union{<:SubOrDArray, <:AbstractDict}, },
																				 <:Tuple{Vararg{<:AbstractDict}}}

													 )::Bool

	for i=2:length(arrays)
		same_mesh_distrib(arrays[1],arrays[i]) || return false 
	end  

	return true 

end 


function findall_indexin(a::AbstractVector, b::AbstractVector
													 )::AbstractVector{Int} 

	(isempty(a)|isempty(b)) && return Int[]

	I = indexin(a,b) 

	i1 = findfirst(!isnothing, I) 

	isnothing(i1) && return Int[]

	i2 = findlast(!isnothing, I)


	isnothing(i2) && return Int[]


	return all(!isnothing, view(I,i1+1:i2-1)) ? (i1:i2) : findall(!isnothing, I)

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



	SignalProcessing.init_run_ms!(storage) 

	setindex!(X, Utils.closest_periodic_shifted_a(X[1], quantized_values, 1), 1)

	SignalProcessing.run_m!(storage, X[1])
	
	for i = 2:lastindex(X)

		setindex!(X, Utils.closest_periodic_shifted_a(X[i], X[i-1], 1), i)

		SignalProcessing.run_m!(storage, X[i])

	end 


	setindex!(S, 
						Utils.closest_periodic_shifted_a(SignalProcessing.run_mean(storage), quantized_values, 1), 
						1)

	setindex!(S, 
						Statistics.stdm(X, Utils.closest_periodic_b(S[1], quantized_values, 1)),
						2)


	return 

end 

function wcc_stat(X::AbstractArray{<:Real,N},
									 quantized_values::AbstractVector{<:Real}...;
									 kwargs...
								 )::Array{Float64,N} where N

	wcc_stat!(copy(X), quantized_values...; kwargs...)

end 


	
function wcc_stat!(X::AbstractArray{<:Real,N},
									 quantized_values::AbstractVector{<:Real}=[0];
									dim::Int=1 
								 )::Array{Float64,N} where N

	# take vectors along axis 'dim' and iterate the other axes 

	storage = SignalProcessing.init_run_m()  # same for all sets of WCCs 


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


function wcc_stat_axpy(a::Real,
											 X::AbstractArray{<:Real},
											 Y::AbstractArray{<:Real},
											 args...;
											 kwargs...)::Vector{Float64}

	wcc_stat_axpy!(a, X, copy(view(Y,:)), args...; kwargs...)

end 		

function wcc_stat_axpy!(a::Real,
											 X::AbstractArray{<:Real},
											 Y::AbstractArray{<:Real},
											 args...;
											 kwargs...)::Vector{Float64}

	wcc_stat!(LinearAlgebra.axpy!(a,view(X,:),view(Y,:)), args...; kwargs...)

end 		

function wcc_stat_diff(args...; kwargs...)::Vector{Float64}

	wcc_stat_axpy(-1, args...; kwargs...)

end 		

function wcc_stat_diff!(args...; kwargs...)::Vector{Float64}

	wcc_stat_axpy!(-1, args...; kwargs...)

end 		

function wcc_stat_sum!(args...; kwargs...)::Vector{Float64}

	wcc_stat_axpy!(1, args...; kwargs...)

end 		

function wcc_stat_sum(args...; kwargs...)::Vector{Float64}

	wcc_stat_axpy(1, args...; kwargs...)

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

#
#function occupied_subspace(k::AbstractVector, H::Function, data...;
#													 kwargs...)#::Tuple{Matrix,Vector}
#
#	occupied_subspace(H(k, data...); kwargs...)
#
#end  
function eigH(k::AbstractVector, H::Function, data...;
													 kwargs...)#::Tuple{Matrix,Vector}

	eigH(H(k, data...); kwargs...)

end  



#function init_occup(args...)
#
#	init_eigH(args...; halfspace=true)
#
#end 

function init_eigH(k::AbstractVector, H::Function, data...; kwargs...) 

	init_eigH(H(k, data...); kwargs...)

end 

function init_eigH(ij::NTuple{2,Int},
									 H::Function, kij::Function, 
									 perturb::AbstractArray{ComplexF64,4}, Hdata...; kwargs...)

	init_eigH(select_mesh_point(perturb,ij); kwargs...)

end 



function init_eigH(ij_to_arg::Function, data...; kwargs...) 

	init_eigH(ij_to_arg(1,1), data...; kwargs...)

end 



function init_eigH(H::AbstractMatrix{T};
									 halfspace::Bool=false,
									 kwargs...
									 )::Vector  where T<:Number 

	n = LinearAlgebra.checksquare(H)

	m = div(n, halfspace+1)

	return [Matrix{promote_type(T,Float64)}(undef, n, m), 
					Vector{Float64}(undef, m)]

end 


#function occupied_subspace(M::AbstractMatrix{T}, args...; kwargs...
#														) where T<:Number 
##														)::Matrix{T} where T<:Number 
#
#	u = init_occup(M)
#
#	occupied_subspace!(u, M, args...; kwargs...)
#
#	return u 
#
#end   

function eigH(M::AbstractMatrix{T}, args...; kwargs...
														) where T<:Number 
#														)::Matrix{T} where T<:Number 

	u = init_eigH(M)

	eigH!(u, M, args...; kwargs...)

	return u 

end  


#function unoccupied_subspace(args...; occupied::Bool=false, kwargs...) 
#
#	occupied_subspace(args...; occupied=false, kwargs...)
#
#end 
#
#function unoccupied_subspace!(args...; occupied::Bool=false, kwargs...) 
#
#	occupied_subspace!(args...; occupied=false, kwargs...)
#
#end 

#function unoccupied_subspace(M::AbstractMatrix{T}, args...; kwargs...
#														) where T<:Number 
##														)::Matrix{T} where T<:Number 
#
#	u = init_occup(M)
#
#	unoccupied_subspace!(u, M, args...; kwargs...)
#
#	return u 
#
#end   

#function eigH!(storage, ij::NTuple{2,Int}, get_k::Function, 
#							 perturb::AbstractArray{ComplexF64,4},
#							 Hdata...; kwargs...)::Nothing 
#
#	eigH!(storage, get_k(ij), Hdata..., select_mesh_point(perturb, ij);
#				kwargs...) 
#
#
#end   

function eigH!(storage::AbstractArray, 
							 ij_or_k::Union{AbstractVector{<:Real}, 
															Tuple{Vararg{<:Real}},
															Int},
							 H::Function,
							 #ij::NTuple{2,Int}, get_k::Function, perturb::AbstractArray{ComplexF64,4},
							 Hdata...; kwargs...)::Nothing 
#	eigH!(d_ij,(i,j),...)

	eigH!!(storage, H(ij_or_k, Hdata...); kwargs...) 


end   
	#

#tuple, H, get_k 


#function eigH!(storage, k::AbstractVector{<:Real}, H::Function, Hdata...;
#														kwargs...)::Nothing 
#
#	eigH!!(storage, H(k, Hdata...); kwargs...)
#
#end  



function eigH!(storage, H::AbstractMatrix; kwargs...)::Nothing 

	eig = LinearAlgebra.eigen(H)

	eigH!(storage, eig.vectors, eig.values; kwargs...)

end  

function eigH!!(storage::AbstractArray, H::AbstractMatrix; kwargs...
							 )::Nothing 


	eig = LinearAlgebra.eigen!(H)

	eigH!(storage, eig.vectors, eig.values; kwargs...)

end  



function eigH!(uv::AbstractVector{<:AbstractArray},
														psi::AbstractMatrix,
														E::AbstractVector{<:Real};
														kwargs...
														)::Nothing 

	copy!.(uv, psien_sorted_energy(psi, E; kwargs...))

	return 

end   


function sortperm_energy(N::Int;
												 halfspace::Bool=false,
												 occupied::Bool=false,
												 kwargs...
												 )::AbstractVector{Int}

	halfspace || return 1:N 

	return occupied ? (1:div(N,2)) : (div(N,2)+1:N)

end 

function sortperm_energy(E::AbstractVector{<:Real}; 
												 halfspace::Bool=false,
												 kwargs...
												 )::AbstractVector{Int}

	halfspace || return sortperm(E; kwargs...)

	return partialsortperm(E, 
												 sortperm_energy(length(E); 
																				 halfspace=halfspace, kwargs...)
												 )
end 

function psien_sorted_energy(
													 psi::AbstractArray{ComplexF64,N},
														E::AbstractVector{<:Real};
														kwargs...
														)::AbstractVector{<:AbstractArray} where N

	@assert N>=2  

	i = sortperm_energy(E; kwargs...)

	return [selectdim(psi, 2, i), view(E, i)]

end   


function psi_sorted_energy(
													 psi::AbstractArray{ComplexF64,N},
													 E::Union{<:Int,<:AbstractVector{<:Real}
																		}=size(psi,2);
												 kwargs...
											 )::AbstractArray{ComplexF64,N} where N

	@assert N>=2 

	return selectdim(psi, 2, sortperm_energy(E; kwargs...))

end 


function eigH!(u::AbstractMatrix{ComplexF64},
							 psi::AbstractMatrix{<:Number},
														E::AbstractVector{<:Real};
														kwargs...
														)::Nothing 

	copy!(u, psi_sorted_energy(psi, E; kwargs...))

	return 

end  



#function occupied_subspace!(uv::AbstractVector{<:AbstractArray},
#														psi::AbstractMatrix,
#														E::AbstractVector{<:Real};
#														kwargs...
#														)::Nothing 
#
#	u,v = uv  
#
#	i = partialsortperm(E, 1:div(length(E),2); kwargs...)
#
#	setindex!(u, selectdim(psi, 2, i), :, :)
#
#	setindex!(v, view(E, i), :)
#
#
#	return 
#
#end 
#
#function occupied_subspace!(
#														u::AbstractMatrix,
#														psi::AbstractMatrix,
#														E::AbstractVector{<:Real};
#														kwargs...
#														)::Nothing 
#
#	occupied_subspace!([u,zeros(size(u,2))], psi, E; kwargs...)
#
#end 
#





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function unitary_part(G::AbstractMatrix)::Matrix

	unitary_part!(copy(G))

end  

function unitary_part!(G::AbstractMatrix)::AbstractMatrix

	svd = LinearAlgebra.svd!(G)  

	LinearAlgebra.mul!(G, svd.U, svd.Vt)

end 

function unitary_overlap(wf1::AbstractMatrix{T1},
								 wf2::AbstractMatrix{T2}
								 )::Matrix{promote_type(T1,T2)
													 } where {T1<:Number,T2<:Number} 

	A = overlap(wf1,wf2)

	unitary_part!(A)

end 

function unitary_overlap!(A::AbstractMatrix{T},
									wf1::AbstractMatrix{T1},
								 wf2::AbstractMatrix{T2}
								 )::AbstractMatrix{T} where {T<:Number,T1<:T,T2<:T}

	overlap!(A, wf1,wf2)

	unitary_part!(A)  

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function check_kMirror_possible(n::Int, k0::Real)::Int 

	@assert n>1 

	m0 = k0*(n-1)/pi 

	m1 = Int(round(m0))

	if !isapprox(m0,m1,atol=1e-10)

		error("Not all k-s will have mirror images in the mesh")

	end  

	return mod1(m1,n-1)

end 

function uniqueInds_kMirror(n::Int, k0::Real)::Vector{Int}

	m = check_kMirror_possible(n, k0)

#m even => origin contained 
#m odd => origin not contained 

#	start = div(m,2) + isodd(m)
#out = mod1.(start+1:start+nr, n-1)

#	@show start=div(m,2)+1+isodd(m)
#	@show div(n,2)+iseven(m)*isodd(n)
#
	out = range(start=div(m,2)+1+isodd(m),length=div(n,2)+iseven(m)*isodd(n)) 

	#@show mod1.(out,n-1) 
#println(get_kij(n,k0)(mod1.(out,n-1)))
	#println(get_kij(n,k0;restricted=false)(out)) 

#	return sort!(mod1.(out,n-1))
	return mod1.(out,n-1)

end 

function uniqueKs_kMirror(n::Int, k0::Real)::Vector{Float64}

	get_kij(n,k0)(uniqueInds_kMirror(n, k0))

end 

function check_nr_kPoints(n::Int)

	@assert n>2 "Too few k points"  

	return n 

end 


function next_kPoint(k_ref::Real, dk::Real 
										 )::Float64

	@assert dk >=0 

	k_ref - dk 

end 


function get_kij!(n::Int, k0::Real; 
								 restricted::Bool=true
								)::Function 

	dk = 2pi/(n-1)
	kmax = 2pi+k0+dk 

	function get_kij_!(k::AbstractVector{Float64},
										 inds::Union{Tuple{Vararg{Int}},
																			<:AbstractVector{Int}}
													 )::AbstractVector{Float64} 

		@assert length(k)==length(inds) 

		for (i,ind)=enumerate(inds)

			restricted && @assert 1<=ind<n

			setindex!(k, next_kPoint(kmax, dk*ind), i)

		end 

		return k 

	end  

	function get_kij_!(k::AbstractVector, inds::Int...
										 )::AbstractVector{Float64}
		get_kij_(inds)

	end 

	return get_kij_!

end    

		
function get_kij(n::Int, k0::Real; 
								 kwargs...
								)::Function 

	kij! = get_kij!(n,k0; kwargs...)

	function get_kij_(inds::Union{Tuple{Vararg{Int}},
																			<:AbstractVector{Int}}
													 )::Vector{Float64} 

		kij!(Vector{Float64}(undef, length(inds)), inds)

	end  

	get_kij_(inds::Int...)::Vector{Float64} = get_kij_(inds)

	return get_kij_ 

end    


function get_kij(n::Int, k0::Real,
								 precomp_ks::AbstractVector{<:Real},precomp_dir::Int; 
								 kwargs...
								)::Function 

	kij = get_kij(n,k0; kwargs...)

	function get_kij_(inds::Union{Tuple{Vararg{Int}},
																			<:AbstractVector{Int}}
													 )::Vector{Float64} 

		k = kij(inds) 

		i,m = fldmod1(inds[precomp_dir], n-1)

#		main period : i=1
#		ind smaller than 0: i=0, -1, ... # nr periods outside 
#    2pi*(1-i) = 2pi, 4pi, etc
# 	ind >= n: i=2, 3, .... # nr periods outside 
#  2pi * (1-i) = -2pi, -4pi, ...
	
		k[precomp_dir] = next_kPoint(precomp_ks[m], 2pi*(i-1))

		return k

	end 

	get_kij_(inds::Int...)::Vector{Float64} = get_kij_(inds)

	return get_kij_ 

end    
function get_kij(n::Int, 
								 precomp_ks::AbstractVector{<:AbstractVector{<:Real}};
								 kwargs...
								)::Function 

	kij = get_kij(n,0; kwargs...)

	function get_kij_(inds::Union{Tuple{Vararg{Int}},
																			<:AbstractVector{Int}}
													 )::Vector{Float64} 

		k = kij(inds) 

		for (precomp_dir,ks) in enumerate(precomp_ks)

			i,m = fldmod1(inds[precomp_dir], n-1)

			k[precomp_dir] = next_kPoint(ks[m], 2pi*(i-1))

		end 

		return k

	end 

	get_kij_(inds::Int...)::Vector{Float64} = get_kij_(inds)

	return get_kij_ 

end    


function get_kij_constrained(n::Int, k0::Real,
									dim::Int, fixed_k::Float64, fixed_dir::Int;
									kwargs...
								)::Function 

	kij! = get_kij!(n,k0; kwargs...)

	inds_k = setdiff(1:2, fixed_dir) 


	function get_kij_line_(inds::Union{Tuple{Vararg{Int}},
																			<:AbstractVector{Int}}
													 )::Vector{Float64} 

		k = Vector{Float64}(undef, dim)

		kij!(view(k,inds_k), inds) 

		setindex!(k, fixed_k, fixed_dir)


	end  

	get_kij_line_(inds::Int...)::Vector{Float64} = get_kij_line_(inds)

	return get_kij_line_ 


end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




#function get_ik(n::Int, k0::Real)::Function 
#
#	check_nr_kPoints(n)
#
#	check_kMirror_possible(n, k0)
#
#
#	return function get_ik_(k::Real)::Int 
#	
#		i = (k0-k)/2pi*(n-1) + 1 
#
#		@assert Int(round(i))≈i
#	
#		return reduce_index(Int(round(i)),n-1)
#
#	end  
#
#end  


function ind_minusk(i::Int, n::Int, k0::Real;
									 restricted::Bool=true)::Int 

	j = Int(round(k0*(n-1)/pi))+2-i

	return restricted ? mod1(j, n-1) : j

end 

function ind_minusk(n::Int, k0::Real;
									 restricted::Bool=true)::Function 

	i0 = check_kMirror_possible(n,k0) + 2 

	restricted && return ind_minusk_r(i::Int)::Int = mod1(i0-i, n-1) 

	return ind_minusk_u(i::Int)::Int = i0-i 

end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function select_mesh_point(dir::Int,
													 data::AbstractArray{T,M},
													 j::Int,
													 )::AbstractArray{T,M-1} where {T<:Number,M}
	
#	selectdim(data,M,j)

	selectdim(data, get_mesh_dir(data, dir), j) 

end 

function select_mesh_point(data::AbstractArray{T,M},
													 j::Int,
													 )::AbstractArray{T,M-1} where {T<:Number,M}
	
	selectdim(data,M,j)

end



function select_mesh_point(data::AbstractArray{T,M},
													 k1::Int,
													 dir1::Int,
													 k2::Int,
													 )::AbstractArray{T,M-2} where {T<:Number,M}

	@assert 1<=dir1<=2 

	select_mesh_point(data,orderinds(dir1,k1,k2))

end   

function select_mesh_point(data::AbstractArray{T,M},
													 i::Union{Int,<:AbstractVector{Int}},
													 j::Union{Int,<:AbstractVector{Int}},
													 )::AbstractArray{T} where {T<:Number,M}

	selectdim(selectdim(data,M,j),M-1,i)

end  

function select_mesh_point(data::AbstractArray{T,M},
													 ij::NTuple{N,Int}
													 )::AbstractArray{T,M-N} where {T<:Number,M,N}

	select_mesh_point(data, ij...)

end 

function select_mesh_point(get_data::Function, 
													 i_or_ij::Union{Int, NTuple{2,Int}},
													 j_or_none::Int...
													 )::AbstractArray

	get_data(i_or_ij, j_or_none...)

end 


function select_mesh_point(
													 data::AbstractVector{<:AbstractArray}, args...
													 )::Vector{<:AbstractArray}
	
	[select_mesh_point(d,args...) for d in data]

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function get_one_wrapper(get_one::Function,
												 data::Union{<:AbstractArray{<:Number},
																		 <:AbstractArray{<:AbstractArray},
																		 },
												 i::Int, 
												 j::Int,
												 args...;
												 kwargs...
												)

	get_one(select_mesh_point(data, i, j), args...; kwargs...)

end 


function get_one_wrapper(get_one::Function,
												 ij_to_arg::Function,
												 i::Int, 
												 j::Int,
												 args...; kwargs...
												)
	
	get_one(ij_to_arg(i,j), args...; kwargs...)

end 

function get_one_wrapper(get_one::Function,
												 data::Union{<:AbstractArray{<:Number},
																		 <:AbstractArray{<:AbstractArray},
																		 },
												 i::Int, 
												 arg1::Union{<:AbstractArray,Float64,Function},
												 args...;
												 kwargs...
												)

	get_one(select_mesh_point(data, i), arg1, args...; kwargs...)

end 


function get_one_wrapper(get_one::Function,
												 i_to_arg::Function,
												 i::Int, 
												 arg1::Union{<:AbstractArray,Float64,Function},
												 args...; kwargs...
												)
	
	get_one(i_to_arg(i), arg1, args...; kwargs...)

end 







function get_one_wrapper!(storage::AbstractArray,
													get_one!::Function,
													ij_to_arg::Function,
													i::Int,
													j::Int,
													data...; 
													kwargs...)

	get_one!(storage, ij_to_arg(i,j), data...; kwargs...)

end 


function get_one_wrapper!(storage::AbstractArray,
													get_one!::Function,
												data::Union{
																		<:SubOrArray{<:Number},
																		<:AbstractVector{<:SubOrArray},
																		<:SubOrSArray{<:Number},
																		<:AbstractVector{<:SubOrSArray},
																		},
												 i::Int, 
												 j::Int,
												 data_...;
												 kwargs...
												)

	get_one!(storage, select_mesh_point(data, i, j), data_...; kwargs...)

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function check_parallel_possible(;parallel::Bool=false, kwargs...)::Bool

	parallel && nworkers()>1 

end  





# internal 
function _inds_distrib_array(array_size::NTuple{N,Int},
														array_workers::AbstractVector{Int},
														array_distrib::AbstractVector{Int};
														kwargs...
														)::Dict{Int, NTuple{N,UnitRange{Int}}
																		 } where N

	
	parallel = check_parallel_possible(;kwargs...)

	I = findall(>(1), array_distrib)

	@assert length(I)==parallel 

	t0 = Tuple(1:a for a=array_size)

	parallel || return Dict(only(array_workers)=>t0)

	i = only(I)

	distr_i = Utils.EqualDistributeBallsToBoxes_cumulRanges(
																array_size[i], array_distrib[i])
	
	return Dict(w=>Tuple(i==d ? J : a for (d,a)=enumerate(t0)) for (w,J)=zip(array_workers,distr_i))

end 
														

function _inds_distrib_array(
														args_as_inst...; kwargs_as_inst...
														)::Dict{Int,Tuple{Vararg{UnitRange{Int}}}}

	_inds_distrib_array(_prep_init_array(args_as_inst...; kwargs_as_inst...)...;
											kwargs_as_inst...)

end 


function inds_distrib_array(args...; 
														custom_ak::Function=init_storage_ak,
														kwargs...
														)::Dict{Int,Tuple{Vararg{UnitRange{Int}}}}

	wrapper_ak(custom_ak, _inds_distrib_array, args...; kwargs...) 

end 




function init_storage(item1::AbstractVector{<:AbstractArray},
											args...; kwargs...)::Vector

	[init_storage(it, args...; kwargs...) for it in item1]

end 


function init_storage(args...; kwargs...)

	wrapper_ak(init_storage_ak, _init_storage, args...; kwargs...)

end 

function init_storage_ida(item1::AbstractVector{<:AbstractArray},
													args...; kwargs...)

	AKs = [init_storage_ak(it, args...; kwargs...) for it=item1]

	return ([_init_storage(a...; k...) for (a,k)=AKs], 
					[_inds_distrib_array(a...; k...) for (a,k)=AKs]
					) 


#	collect(zip((init_storage_ida(it, args...; kwargs...) for it in item1)...))

end   

function init_storage_ida(args...; 
													custom_ak::Function=init_storage_ak,
													kwargs...)

	wrapper_ak(custom_ak, (_init_storage, _inds_distrib_array), args...; kwargs...)

end   




#===========================================================================#
#
# _prep_init_array and _init_storage take identical same args and kwargs
# internal methods 
#
#---------------------------------------------------------------------------#


function _prep_init_array(
											s::NTuple{Ns,Int},
											ns::NTuple{Nn,Int},
											args...;
											distr_axis::Int=Nn,
											kwargs...
											)::Tuple{NTuple{Ns+Nn,Int}, Vector{Int}, Vector{Int}
															 } where {Ns,Nn}
	
	array_size = (s..., (n-1 for n=ns)...)

	d = ones(Int, Ns+Nn)


	parallel = check_parallel_possible(;kwargs...)



	parallel || return (array_size, [myid()], d) # just the master proc

	@assert 1<=distr_axis<=Nn
	
	w = workers()    

	wmax = min(length(w), array_size[Ns+distr_axis])

	return (array_size, w[1:wmax], setindex!(d, wmax, Ns+distr_axis))

end 

function _init_storage(#T::DataType, 
											s::NTuple{Ns,Int},
											ns::NTuple{Nn,Int},
											T::DataType;
											parallel::Bool=false,
											shared::Bool=false,
											kwargs...
											)::AbstractArray{<:Number,Ns+Nn} where Ns where Nn 

	array_size, array_workers, array_distrib = _prep_init_array(s,ns; 
																						parallel=parallel, 
																						kwargs...)


	parallel || return zeros(T, array_size)

	shared || return dzeros(T, array_size, array_workers, array_distrib)



#	p = mkpath(abspath("Data/SharedArrays"))

#	fn = joinpath(p,replace(string(mod(time(),1e6)),"."=>"")*".tmp")

	return SharedArray{T}(#fn,
												array_size; 
												pids=union(myid(), array_workers),
#												mode="w+",
												)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#







function init_storage_ak(T::DataType, ns::Int...; kwargs...)

	init_storage_ak(T, (), ns...)

end 

function init_storage_ak(T::DataType, 
												 s::Tuple{Vararg{Int}}, 
												 n::Vararg{Int}; 
												 kwargs...)

	init_storage_ak(T, s, n; kwargs...)

end 

function init_storage_ak(T::DataType, 
												 s::Tuple{Vararg{Int}}, 
												 n::Tuple{Vararg{Int}}; 
												 kwargs...)

	(s, n, T), kwargs

end 

function init_storage_ak(
												 s::Tuple{Vararg{Int}}, 
												 n::Tuple{Vararg{Int}},
												 T::DataType;
												 kwargs...)

	(s, n, T), kwargs

end 




function init_storage_ak(item1::AbstractArray{T}, n::Int; kwargs...) where T

	init_storage_ak(T, size(item1), n, n; kwargs...)

end 

function init_storage_ak(item1::AbstractArray{T}, 
												 ns::Tuple{Vararg{Int}}; kwargs...) where T

	init_storage_ak(T, size(item1), ns; kwargs...)

end 




function ak_from_storage(size_array_mesh::NTuple{N,Int},
												 T::DataType...;
												 kwargs...
												) where N

	@assert N>=2 

	(size_array_mesh[1:2], 
	 Tuple(nr_kPoints_from_mesh(size_array_mesh,d,N-2) for d=1:N-2),
	 T...
	 ), kwargs


end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function store_on_mesh1(get_one::Function, get_one!::Function,
												source::Union{<:AbstractArray{<:Number},
																			<:AbstractVector{<:AbstractArray},
																		},
												data...; kwargs...
												)::AbstractArray 

	store_on_mesh1(get_one, get_one!, nr_kPoints_from_mesh1(source),
								source, data...; kwargs...)

end 


function store_on_mesh1(get_one::Function, get_one!::Function,
											 n::Int,
												source,
											 data...
												)

#or :	source(1)

	dest = init_storage(get_one_wrapper(get_one, source, 1, data...), (n,))

	store_on_mesh1!!(get_one!, n, source, dest, data...)

end   



#function store_on_mesh(get_one::Function, 
#											 n::Int,
#												source::Union{Function, 
#																		AbstractArray{<:Number},
#																		NTuple{N, <:AbstractArray} where N,
#																		},
#											 data...;
#											 kwargs...
#												)
#
#	error() 
#
#	dest = init_storage(get_one_wrapper(get_one, source, 1,1, data...; 
#																				 kwargs...), n) 
#
#
#	store_on_mesh!(get_one, n, source, dest, data...; kwargs...)
#
#	return dest
#
#end  
function store_on_mesh(get_one::Function, get_one!::Function,
											 n::Int,
											 source,
											 data...;
											 kwargs...
											)

	dest, array_distrib = init_storage_ida(
																	get_one_wrapper(get_one, source, 1,1, data...), 
																	n; 
																	kwargs...)

	store_on_mesh!!(get_one!, n, source, dest, data...; 
									array_distrib=array_distrib)


end  

#function store_on_mesh(get_one::Function, get_one!,
#												source::AbstractArray{<:Number},
#											 data...; kwargs...
#												)
#	error()
#
#	store_on_mesh(get_one, nr_kPoints_from_mesh(source), source, data...;
#								kwargs...)
#
#end   

function store_on_mesh(get_one::Function, get_one!::Function,
												source::Union{<:AbstractArray{<:Number},
																			<:AbstractVector{<:AbstractArray},
																		},
												data...; kwargs...
												)::AbstractArray 

	store_on_mesh(get_one, get_one!, nr_kPoints_from_mesh(source),
								source, data...; kwargs...)

end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function store_on_mesh!!(get_one!::Function, 
												source::Function, 
												dest::AbstractArray,
												data...; kwargs...
												)::AbstractArray

	store_on_mesh!!(get_one!, nr_kPoints_from_mesh(dest),
									source, dest, data...; kwargs...) 

end 

function store_on_mesh!!(get_one!::Function, 
												source::AbstractArray,
												dest::AbstractArray=source,
												data...; kwargs...
												)::AbstractArray

	store_on_mesh!!(get_one!, nr_kPoints_from_mesh(source),
									source, dest, data...; kwargs...)

end 



function store_on_mesh!!(get_one!::Function, n::Int, args...; 
												 kwargs...)::AbstractArray

	store_on_mesh!!(get_one!, (1:n-1,1:n-1), args...; kwargs...)

end 

function _store_on_mesh!!(get_one!::Function, 
												 inds::Tuple{<:AbstractVector{Int},
																		 <:AbstractVector{Int}},
												source::Union{<:Function, 
																			<:SubOrArray{<:Number},
																			<:AbstractVector{<:SubOrArray},
																			<:SubOrSArray{<:Number},
																			<:AbstractVector{<:SubOrSArray},
																		},
												dest::Union{<:SubOrArray{<:Number},
																		<:SubOrSArray{<:Number},
																		<:AbstractVector{<:AbstractArray},
																		},
												data...; 
												parallel::Bool=false,
												array_distrib=nothing,
												kwargs...
												)::AbstractArray

	for j=inds[2], i=inds[1]

		store_on_mesh_one!!(get_one!, source, i, j, dest, data...; kwargs...)

	end 

	return dest 

end    

function store_on_mesh!!(get_one!::Function, 
												 inds::Tuple{<:AbstractVector{Int},
																		 <:AbstractVector{Int}},
												source::Union{<:Function, 
																			<:SubOrArray{<:Number},
																			<:AbstractVector{<:SubOrArray},
																		},
												dest::Union{<:SubOrArray{<:Number},
																		<:AbstractVector{<:SubOrArray},
																		},
												data...; 
												kwargs...
												)::AbstractArray

	_store_on_mesh!!(get_one!, inds, source, dest, data...; kwargs...)

end   






function store_on_mesh!!(get_one!::Function, 
												 inds::Tuple{<:AbstractVector{Int},
																		 <:AbstractVector{Int}},
												source::Union{<:Function, 
																			<:AbstractArray{<:Number},
																			<:AbstractVector{<:AbstractArray},
																		},
												 dest::Union{<:AbstractVector{<:SubOrDArray},
																		 SubOrDArray{<:Number},
																		 },
												 # both source and dest may be overwritten
												data...; # assumes data is read-only 
												array_distrib=nothing,
												kwargs...
												)::AbstractArray

	map(procs(eltype(dest)<:AbstractArray ? dest[1] : dest)) do p

		@spawnat p begin   

			store_on_mesh!!(get_one!, 
											map(findall_indexin, get_big_inds(dest), inds),
											localpart_(source,dest),
											localpart_(dest),
										data...;
										kwargs...)

#eventually: get_one!(dest_ij, src_ij, data...) 

		end  

	end .|> fetch 

	return dest 

end   

function store_on_mesh!!(get_one!::Function, 
												 inds::Tuple{<:AbstractVector{Int},
																		 <:AbstractVector{Int}},
												source::Union{<:Function, 
																		 <:AbstractVector{<:SubOrSArray},
																		 SubOrSArray{<:Number},
																		},
												 dest::Union{
																		 <:AbstractVector{<:SubOrSArray},
																		 SubOrSArray{<:Number},
																		 },
												data...; #
												array_distrib::Dict{Int,Tuple{Vararg{UnitRange{Int}}}},
												kwargs...
												)::AbstractArray

	for d in data 

		@assert !isa(d,SubOrDArray)

		isa(d,SubOrSArray) && continue 
	
		@assert !isa(d, SubOrDArray)

		isa(d,SubOrArray) || continue 
		
		length(d) < max(100, sum(length, dest)/50) && continue 

		@warn "Large non-shared array passed? $(size(d))"

		@show typeof(d)
	end 

	map(procs_(array_distrib)) do p 

		locinds = array_distrib[p]

		@spawnat p begin   

			_store_on_mesh!!(get_one!, 
											map(findall_indexin, get_big_inds(locinds), inds),
											localpart_(source,locinds),
											localpart_(dest,locinds),
										data...;
										kwargs...)

		end  

	end .|> fetch 

	return dest 

end   










#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function store_on_mesh1!!(get_one!::Function, 
												source::Function,
												dest::SubOrArray{<:Number},
												data...; kwargs...
												)::AbstractArray

	store_on_mesh1!!(get_one!, nr_kPoints_from_mesh1(dest), 
									 source, dest, data...; kwargs...)
end 

function store_on_mesh1!!(get_one!::Function, 
												source::SubOrArray{<:Number},
												args...; kwargs...
												)::AbstractArray 

	store_on_mesh1!!(get_one!, nr_kPoints_from_mesh1(source), 
									 source, args...; kwargs...)

end  


function store_on_mesh1!!(get_one!::Function, 
													nk::Int,
												source::Union{<:Function, 
																			<:SubOrArray{<:Number},
																			<:AbstractVector{<:SubOrArray},
#																			<:SubOrSArray{<:Number},
#																			<:AbstractVector{<:SubOrSArray},
																		},
												dest::Union{<:SubOrArray{<:Number},
																		<:AbstractVector{<:SubOrArray}
																		},
												data...; kwargs...
												)::AbstractArray

	for i=1:nk-1

		store_on_mesh1_one!!(get_one!, source, i, dest, data...; kwargs...)

	end 

	return dest 

end  



#function store_on_mesh1!(get_one::Function, 
#												n::Int,
#												source::Union{Function, AbstractArray{<:Number},
#																		NTuple{N, <:AbstractArray} where N,
#																		},
#												dest::AbstractArray,
#												data...; kwargs...
#												)::Nothing 
#
#	for i=1:n-1
#
#		store_on_mesh1_one!(get_one, source, i, dest, data...; kwargs...)
#
#	end 
#
#	return 
#
#end 
#


function store_on_mesh1_one!!(
														 get_one!::Function, 
												i_to_arg::Function,
														i::Int, 
														dest::Union{<:AbstractArray{<:Number,N},
																					 <:AbstractVector{<:AbstractArray}
																					 },
														
												data...; kwargs...
												)::Nothing where N
	
	get_one!(select_mesh_point(dest,i), i_to_arg(i), data...; kwargs...) 


	return 

end  

function store_on_mesh1_one!!(
														 get_one!::Function, 
														 source::AbstractArray{<:Number},
														i::Int, 
														dest::Union{<:AbstractArray{<:Number,N},
																					 <:AbstractVector{<:AbstractArray}
																					 }=source,
														
												data...; kwargs...
												)::Nothing where N
	
	get_one!(select_mesh_point(dest,i), select_mesh_point(source,i), data...; 
					 kwargs...) 

	return 

end  






function store_on_mesh_one!!(
							get_one!::Function, 
							source,
							i::Int, j::Int,
							dest::Union{<:SubOrArray{<:Number},
#													<:AbstractVector{<:SubOrArray},
													<:SubOrSArray{<:Number},
#													<:AbstractVector{<:SubOrSArray},
													<:AbstractVector{<:AbstractArray},
													},

							data...; kwargs...
												)::Nothing 


	get_one_wrapper!(select_mesh_point(dest, i, j),
									 get_one!, 
									 source, i, j, 
									 data...; kwargs...)

	return 

end  




function store_on_mesh_one!(get_one::Function, 
												source::Union{Function, AbstractArray{<:Number},
																		NTuple{M, <:AbstractArray} where M,
																		},
														i::Int, j::Int,
												dest::AbstractArray{<:Number},
												data...; kwargs...
												)::Nothing 
	
	copy!(select_mesh_point(dest, i, j),
				get_one_wrapper(get_one, source, i, j, data...; kwargs...)
				)

	return 

end  
function store_on_mesh1_one!(get_one::Function, 
												source::Function,
														i::Int, 
												dest::AbstractArray{<:Number},
												data...; kwargs...
												)::Nothing 
	
	copy!(select_mesh_point(dest, i),
				get_one(source(i), data...; kwargs...)
				)

	return 

end  



function store_on_mesh_one!(get_one::Function, 
												source::Union{Function, AbstractArray{<:Number}},
														i::Int, j::Int,
														dest::AbstractVector{<:AbstractArray},
												data...; kwargs...
												)::Nothing 

	for (k,item) in enumerate(get_one_wrapper(get_one, source, i, j, data...;
																					 kwargs...))

		copy!(select_mesh_point(dest[k], i, j), item)

	end 

	return 

end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_item_ij(ij::Union{#AbstractVector{Int}, 
															 Int, NTuple{2,Int}},
										 get_item_::Function,
										 )::AbstractArray

	get_item_(ij) 

end 

function get_item_ij(#(i,j)::Union{AbstractVector{Int}, Tuple{Int,Int}},
										 ij::NTuple{M,Int},
										 storage::AbstractArray{T,N}
										 )::AbstractArray{T,N-M} where {M,N,T<:Number}

	get_item_ij(ij..., storage)

end 
function get_item_ij(j::Int,
										 storage::AbstractArray{T,N}
										 )::AbstractArray{T,N-1} where T<:Number where N

	PeriodicFuns.cycseldim(storage, N, j)

end  

function get_item_ij(i::Int, j::Int, 
										 storage::AbstractArray{T,N}
										 )::AbstractArray{T,N-2} where T<:Number where N

	PeriodicFuns.cycseldim(PeriodicFuns.cycseldim(storage, N, j), N-1, i)

end  

function get_item_ij(storage::AbstractArray{T}
										 )::Function where T<:Number 

	get_item_ij_(args...)::AbstractArray{T} = get_item_ij(args..., storage)

end 


function get_item_ij(k::Int, dir::Int,
#										 start_ij::Union{AbstractVector{Int},NTuple{2,Int}},
										 start::NTuple{M,Int},
										 storage::AbstractArray{T,N}
										 )::AbstractArray{T,N-M} where {T<:Number,N,M}


	get_item_ij(Tuple(d==dir ? s+k-1 : s for (d,s)=enumerate(start)), storage)

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

function wlo1_from_mesh1(Bloch_WFs::AbstractArray{ComplexF64,3},
											 )::Matrix{ComplexF64}

	mapreduce(*,axes(Bloch_WFs,3)) do n 

		unitary_overlap(get_item_ij(n, Bloch_WFs), get_item_ij(n+1, Bloch_WFs),)

	end 

end  

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

#
#function matmulpairs!(dst::AbstractArray{T,3},
#											src::AbstractArray{T,3},
#											N::Int=size(src,3)
#											)::AbstractMatrix{T} where T<:Number 
#
#	N == 1 && return selectdim(src,3,1)
#	
#	offset = isodd(N) 
#
#	offset && copy!(selectdim(dst,3,1), selectdim(src,3,1))
#
#	for (i,i2) in zip(1+offset:div(N,2)+offset,2+offset:2:N)
#
#		LinearAlgebra.mul!(selectdim(dst,3,i),
#											 selectdim(src,3,i2-1),
#											 selectdim(src,3,i2)
#											 )
#	end 
#
#	matmulpairs!(src, dst, div(N+1,2))
#
#end 

function matmulpairs!(aux::AbstractArray{T,3},
											matrices::Vararg{<:AbstractMatrix{T},N}
											)::AbstractMatrix{T} where T<:Number where N 

	for i=1:div(N,2) 

		LinearAlgebra.mul!(selectdim(aux,3,i),matrices[2i-1],matrices[2i])

	end 

	isodd(N) && copy!(selectdim(aux,3,div(N+1,2)), matrices[N])

	return matmulpairs!(aux, div(N+1,2))

end 

function matmulpairs!(dst::AbstractArray{T,3},
											src::AbstractArray{T,3},
											inds_src::AbstractVector{Int},
											i0::Int...
											)::AbstractMatrix{T} where T<:Number  

	matmulpairs!(dst, selectdim(src,3,inds_src), i0...)

end 


function matmulpairs!(dst::AbstractArray{T,3},
											src::AbstractArray{T,3},
											i0::Int=1,
											)::AbstractMatrix{T} where T<:Number 
	
	Npairs,offset = divrem(size(src,3),2)

	Npairs==0 && return selectdim(src,3,1)

	offset==1 && copy!(selectdim(dst,3,i0), selectdim(src,3,1)) 
									# not touched at this iteration 

# dst[i0:imax] will contain the results 
	imax = i0+Npairs+offset-1

	for (i,i2) in zip(i0+offset:imax, 2+offset:2:size(src,3))

		LinearAlgebra.mul!(selectdim(dst,3,i),
											 selectdim(src,3,i2-1),
											 selectdim(src,3,i2)
											 )
	end 

	return matmulpairs!(dst, dst, i0:imax, 1+imax*(size(dst,3)-i0+1>imax))

end 



function matmulpairs!(data::AbstractArray{T,3}, # overwritten 
											 stop::Int 
											)::AbstractMatrix{T} where T<:Number 

#	stop = size(data,3)-1

	N::Int = stop 

	N>1 && @assert stop<size(data,3) "Extra space needed for storage" 


	while N>1 

		for (j,k) in zip(stop:-1:stop-div(N,2)+1, stop-2:-2:stop-N)
	
			LinearAlgebra.mul!(selectdim(data, 3, mod(j,size(data,3))+1),
												 selectdim(data, 3, mod(k,size(data,3))+1),
												 selectdim(data, 3, mod(k+1,size(data,3))+1),
												 )
	
		end 
	
	
		if isodd(N) 
	
			copy!(selectdim(data, 3, mod(stop-div(N,2), size(data,3))+1),
						selectdim(data, 3, mod(stop-N, size(data,3))+1))
	
		end  
	
		stop = mod(stop,size(data,3))+1
	
		N = div(N+1,2)

	end 
	
	
	return selectdim(data,3,stop)


end  


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function unitary_overlaps_line!(
#												overlaps::AbstractArray{ComplexF64,3},
#												dir::Int,
#												i_perp::Int,
#												wfs::AbstractArray{ComplexF64,4},
#												nmesh::Int...
#											 )::AbstractArray{ComplexF64,3}
#
#	unitary_overlaps_line!(overlaps, select_mesh_point(3-dir, wfs, i_perp),
#												 nmesh...)
#
#end 





function unitary_overlaps_line!(
												overlaps::AbstractArray{ComplexF64,3},
												wfs::AbstractArray{ComplexF64,3},
												nmesh::Int=nr_kPoints_from_mesh1(wfs),
											 )::AbstractArray{ComplexF64,3}

	for k = 1:nmesh-1  

		unitary_overlap!(selectdim(overlaps,3,k),
										 get_item_ij(k, wfs),
										 get_item_ij(k+1, wfs),
										 )

	end  

	return overlaps 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



# obsolete 
#function wlo_from_mesh!(dst::AbstractMatrix{ComplexF64},
#												IJ0::Union{AbstractVector{Int},Tuple{Int,Int}},
#												dir::Int,
#											 wfs::AbstractArray{ComplexF64,4},
#											 n::Int,
#											 overlaps::AbstractArray{ComplexF64,3},
#											 )::Nothing 
#
#	copy!(dst, matmulpairs!(unitary_overlaps_line!(overlaps, n, dir, IJ0, wfs)))
#
#	return  
#
#######


#wlo1_one_inplace!(
#																overlaps
#																wfs 
#																dir, IJ0 
#
#	wlo_from_ov!(unitary_overlaps_line!(overlaps, 
#																				 prep_args_uoom(args_uoom...)...,
#																				 WFs), 
#							 ov_aux,
#							 )

#####

#end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function init_arrays_gap(nk::Int, h::Function, Hdata...)::Tuple

	WF = init_storage(init_eigH(rand(2), h, Hdata...)[1], (nk,))

	WF_occ = psi_sorted_energy(WF; halfspace=true, occupied=true) 
						# only to init the right size 

	return (WF, init_overlaps_line(WF_occ))

end 
function pack_data_gap(Hdata, (nk,k0), sector)

		(
							Hdata,
							(nk, k0),
							init_arrays_gap(nk, Hdata...),
							sector,
							)

end 




""" 
Calculate the gap of W_dir1 for parameter k[dir2]=k2
"""
function calc_gap!(
									 (Hdata,
										(nk,k0),
										(WF, (overlaps,ov_aux)), 
										(occupied,dir1),
										),
							k2::Float64
							)::Float64

	get_K = get_kij_constrained(nk, k0, 2, k2, 3-dir1)
	
	store_on_mesh1!!(eigH!, nk, get_K, WF, Hdata...)

	WF_occ = psi_sorted_energy(WF; halfspace=true, occupied=occupied)

	w1 = wlo1_one_inplace!(overlaps, ov_aux, WF_occ)

	return Wannier_min_gap(get_periodic_eigvals!(w1))

end  


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function psiH_on_mesh1(n::Int, k0::Real, H::Function, Hdata...; kwargs...
											)::Array{ComplexF64,3}

	@assert !get(kwargs,:parallel,false) "Not implemented"

	store_on_mesh1(first∘init_eigH, eigH!, n, 
								 get_kij(n,k0), H, Hdata...; kwargs...)

end  

function psiH_on_mesh1(n::Int, k0::Real, 
											perturb::AbstractArray{ComplexF64,3},
											 H::Function, Hdata...; kwargs...
											)::Array{ComplexF64,3}

	@assert !get(kwargs,:parallel,false) "Not implemented" 

	kij = get_kij(n,k0)  

	WFs = init_storage1(init_eigH(kij(1), H, Hdata...; kwargs...)[1], n; 
													 kwargs...)

	store_on_mesh1!!(eigH!, n, identity, WFs, H, kij, perturb, Hdata...; kwargs...)

	return WFs 

end  


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function psiH_on_mesh(n::Int, kij::Function, H::Function, Hdata...; 
											kwargs...
											)::Union{Array{ComplexF64,4},SharedArray{ComplexF64,4}} 

	out = store_on_mesh(first∘init_eigH, eigH!, n, kij, H, Hdata...; kwargs...)

	return out isa SubOrDArray ? convert(Array, out) : out 

end  
function psiH_on_mesh(n::Int, 
											k::Union{<:Real,
															 <:AbstractVector{<:AbstractVector{<:Real}}
															 },
											args...; kwargs...
											)::AbstractArray{ComplexF64,4}

	psiH_on_mesh(n, get_kij(n,k), args...; kwargs...)

end 


function psiH_on_mesh(n::Int, 
											kij::Function,
											perturb::AbstractArray{ComplexF64,4},
											H::Function, Hdata...; kwargs...
											)::Union{Array{ComplexF64,4},SharedArray{ComplexF64,4}} 

	out = store_on_mesh(first∘init_eigH, eigH!, 
								n, tuple, H, kij, perturb, Hdata...; kwargs...) 

	return out isa SubOrDArray ? convert(Array, out) : out 
##
#	store_on_mesh!!(eigH!, n, source=tuple, dest=WFs,  rest...)
#	rest= (H,kij,perturb,Hdata...)
#	Hdata = (bs,pert_strength)
#
#	store_one!!(eigh!, source=tuple, i, j, dest=WFs, rest...)
#
#	get_one_wrapper!(WFs[i,j], eigH!, tuple, i, j, rest...)
#
#	eigH!(WFs[i,j], (i,j), rest...) 
#	
#	eigH!(WFs[i,j], (i,j), H, rest...) 
#	rest= (kij,perturb,Hdata...)
#
# eigH(WFs[i,j], H(ij, rest...); kwargs...)
#
# H(ij, kij, perturb, bs, pert_strength)
#
##
end 






function enpsiH_on_mesh(n::Int, k0::Real, H::Function, Hdata...; kwargs...
												)::Vector{Array}


	store_on_mesh(init_eigH, eigH!, n, get_kij(n,k0), H, Hdata...; kwargs...)

end 

function enpsiH_on_mesh(n::Int, k0::Real, 
											perturb::AbstractArray{ComplexF64,4},
											H::Function, Hdata...; kwargs...
												)::Vector{Array}

	kij = get_kij(n,k0)   

	storage = init_storage(init_eigH(kij, H, Hdata...; kwargs...), n; kwargs...)

	store_on_mesh!!(eigH!, n, tuple, storage, H, kij, perturb, Hdata...; kwargs...) 

	return storage 

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function has_symm_on_mesh_one(ij::NTuple{2,Int},
															repres::NTuple{2,Function}, 
															A::AbstractArray{ComplexF64,4};
															kwargs...
															)::BitVector  

	has_symm_on_mesh_one(ij,  repres, get_item_ij, A; kwargs...)
							
end  

function has_symm_on_mesh_one(ij::NTuple{2,Int},
															(op,fij)::NTuple{2,Function}, 
															# act on H and k spaces, respectively
															get_array::Function,
															data...; kwargs...
															)::BitVector  

	BitVector((Symmetries.has_symm(op, 
																 get_array(ij, data...), 
																get_array(fij(ij), data...);
																kwargs...),
						 ))
							
end  



function has_symm_on_mesh(n::Int, k0::Real, # k0 not used 
													symm::NTuple{2,Function},
													A::AbstractArray;
													kwargs...
													)::BitArray{3}

	error() 

	store_on_mesh(has_symm_on_mesh_one, n, tuple, symm, A; kwargs...)


end 


function has_symm_on_mesh(n::Int, k0::Real, symm::NTuple{2,Function},
													H::Function, args...; kwargs...
													)::BitArray{3}

	error() 

	store_on_mesh(has_symm_on_mesh_one, 
								n, tuple, symm, H, get_kij(n,k0), args...;
								kwargs...)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function symmetrize_on_mesh!(A::AbstractArray{ComplexF64,4},
														 inds::Union{Int,Tuple{<:AbstractVector{Int},
																									 <:AbstractVector{Int}}},
														 (op,fij)::NTuple{2,Function}; # symm;
														 kwargs...
														 )::AbstractArray

	store_on_mesh!!(symmetrize_on_mesh_one!, inds, fij, A, A, op; kwargs...)

end 

function symmetrize_on_mesh_one!(Aij::AbstractMatrix{ComplexF64},
									opij::NTuple{2,Int}, 
									seed::AbstractArray{<:Number,4},
									op::Function; kwargs...)::Nothing

	Symmetries.symmetrize!!(Aij, select_mesh_point(seed, opij), op;
													kwargs...)

end 





function wlo1_on_mesh(dir1::Int, Bloch_WFs::AbstractArray,
											)::Array

	error() 

	store_on_mesh(wlo_from_mesh, 
								nr_kPoints_from_mesh(Bloch_WFs), 
								tuple, 
								dir1, 
								Bloch_WFs) 


end  


function orderinds(dir1::Int, i::Union{Int,NTuple{1,Int}})::NTuple{1,Int}

	@assert dir1==1 

	Tuple(i)

end 

	
function orderinds(dir1::Int, (i,j)::NTuple{2,Int})::NTuple{2,Int}

	orderinds(dir1, i, j)

end 

function orderinds(dir1::Int, i::Int, j::Int)::NTuple{2,Int}

	dir1==1 && return (i,j)

	dir1==2 && return (j,i)

	error()

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#
function loop_inds!(inds::AbstractVector{Int},
										start::Int,
										nmesh::Int=length(inds)+1,
										)::AbstractVector{Int}

	setindex!(inds, start:nmesh-1, 1:nmesh-start)

	setindex!(inds,1:start-1, nmesh-start+1:nmesh-1)

	return inds 

end  

function loop_inds(start::Int, nmesh::Int
										)::AbstractVector{Int}

	loop_inds!(init_storage(Int, nmesh), start, nmesh)

end  




function wlo_from_ov!(
							overlaps::AbstractArray{ComplexF64,3},
							ov_aux::AbstractArray{ComplexF64,3},		 		# mod
							start::Int		
							)::AbstractMatrix{ComplexF64}

	wlo_from_ov!(overlaps, ov_aux, loop_inds(start,nr_kPoints_from_mesh1(overlaps)))

end 
#function wlo_from_ov!(
#							overlaps::AbstractArray{ComplexF64,3},			# mod 
#							start::Int		
#							)::AbstractMatrix{ComplexF64}
#
#	wlo_from_ov!(overlaps, loop_inds(start,nr_kPoints_from_mesh1(overlaps)))
#
#end 

function wlo_from_ov!(
							overlaps::AbstractArray{ComplexF64,3},
							ov_aux::AbstractArray{ComplexF64,3}, 				# mod 
							ov_inds::AbstractVector{Int}...
							)::AbstractMatrix{ComplexF64}

	matmulpairs!(ov_aux, overlaps, ov_inds...)

end  



#function wlo_from_ov!(
#							overlaps::AbstractArray{ComplexF64,3}, 			# mod 
#							ov_inds::AbstractVector{Int}...
#							)::AbstractMatrix{ComplexF64}
#
#	matmulpairs!(overlaps, ov_inds...)
#
#end  




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#







"""Compute overlaps on a line and and multiply them"""
function wlo1_one_inplace!(overlaps::AbstractArray{ComplexF64,3},		# modif
													 ov_aux::AbstractArray{ComplexF64,3},			# modified
													 WFs::AbstractArray{ComplexF64,3},
													)::AbstractMatrix{ComplexF64} 

	wlo_from_ov!(unitary_overlaps_line!(overlaps, WFs), ov_aux)

end 



# time complexity O(N^2)
function wlo1_on_mesh1_inplace!(
																dest::AbstractArray{ComplexF64,4},				#
																overlaps::AbstractArray{ComplexF64,3},		#
																ov_aux::AbstractArray{ComplexF64,3},			#

																WFs::AbstractArray{ComplexF64,4},

																dir::Int, k_perp::Int,

																)::AbstractArray{ComplexF64,4}


	wlo1_on_mesh1_inplace!(select_mesh_point(3-dir, dest, k_perp),
												 overlaps, ov_aux,
												 select_mesh_point(3-dir, WFs, k_perp)
												 )

	return dest 

end 

function wlo1_on_mesh1_inplace!(
																dest::AbstractArray{ComplexF64,3},				#
																overlaps::AbstractArray{ComplexF64,3},		#
																ov_aux::AbstractArray{ComplexF64,3},			#
																WFs::AbstractArray{ComplexF64,3},
																)::AbstractArray{ComplexF64,3}


	copy!(select_mesh_point(dest, 1),
				wlo1_one_inplace!(overlaps, ov_aux, WFs), # computes overlaps first 
				)

	# time complexity O(N^0.9), a bit worse for larger matrices
	for m=2:nr_kPoints_from_mesh1(overlaps)-1

		copy!(select_mesh_point(dest, m),
						matmulpairs!(ov_aux,
											 inv(select_mesh_point(overlaps, m-1)),
											 select_mesh_point(dest, m-1),
											 select_mesh_point(overlaps, m-1),
											 )
					)
									
					#wlo_from_ov!( overlaps, ov_aux, m), # MUCH more expensive !

	end 

	return dest 

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function wlo_filter_distr_kwargs(dir::Int;
																distr_axis::Int=10,
																kwargs...
																)

	(distr_axis=3-dir, kwargs...)

end  

function wlo_filter_distr_kwargs(;
																 dir::Int,
																distr_axis::Int=10,
																kwargs...
																)

	wlo_filter_distr_kwargs(dir; kwargs...)

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function init_overlaps_line_ak(nwf::Int, nmesh::Int; 
											parallel::Bool=false,
											distr_axis::Int=6,
											kwargs...
											)

	parallel && @assert 1<=distr_axis<=2 



	ns = parallel ? Tuple(i==distr_axis ? min(nworkers()+1,nmesh) : nmesh for i=1:2) : (nmesh,)

	return init_storage_ak(ComplexF64, (nwf,nwf), ns;
												 parallel=parallel, distr_axis=distr_axis, kwargs...)

end 

function init_overlaps_line_ak(wfs::AbstractArray{ComplexF64,3},
													 )

	init_overlaps_line_ak(1,wfs)

end  

function init_overlaps_line_ak(
														wfs::AbstractArray{ComplexF64,N},
														nmesh::Int,
													 ) where N

	@assert 3<=N<=4 

	init_overlaps_line_ak(size(wfs,2), check_nr_kPoints(nmesh))

end 




function init_overlaps_line_ak(
														dir::Int,
														s::NTuple{N,Int};
														kwargs...
														) where N

	@assert 3<=N<=4

	(s,n,), = ak_from_storage(s) 

	return init_overlaps_line_ak(dir, s, n; kwargs...) 

end 

function init_overlaps_line_ak(dir::Int,
													s::NTuple{2,Int},
													n::NTuple{D,Int},
													T...
													;
													parallel::Bool=false,
											 kwargs...
											) where D 

	@assert 1<=dir<=D<=2 

#	init_wlo_mesh_ak(s[2], n...; wlo_filter_distr_kwargs(dir; kwargs...)...)

	init_overlaps_line_ak(s[2], n[dir];
												parallel=D==2 && parallel,
												wlo_filter_distr_kwargs(dir; kwargs...)...
												)

end  



function init_overlaps_line_ak(
														dir::Int, wfs::AbstractArray{ComplexF64};
														kwargs...
														)

	init_overlaps_line_ak(dir, size(wfs); kwargs...)

end 

function _init_overlaps_line(a...; k... 
														)::Vector{<:Union{Array,DArray,SharedArray}}

	[init_storage(a...; k...), init_storage(a...; k...)]

#	wrapper_ak(2, init_overlaps_line_ak, init_storage, args...; kwargs...)
	
end 

function init_overlaps_line(args...; kwargs... 
														)::Vector{<:Union{Array,DArray,SharedArray}}

	wrapper_ak(init_overlaps_line_ak, _init_overlaps_line, args...; kwargs...)
	
end  

function init_overlaps_line_ida(args...; kwargs... 
														)#::Vector{<:Union{Array,DArray,SharedArray}}

	wrapper_ak(init_overlaps_line_ak, 
						 (_init_overlaps_line, _inds_distrib_array), args...; kwargs...)

end 

function Wannier_sector_storage(W::AbstractArray{T,N}
															)::AbstractArray{T,N} where {T<:Number,N}

	nr_wf = size(W,1)
	
	@assert size(W,2)==nr_wf 

	I = 1:div(nr_wf,2)
	
	return selectdim(selectdim(W,1,I),2,I)

end 



#function init_overlaps_line(
#														dir::Int,
#														wfs::AbstractArray{ComplexF64,N};
#														parallel::Bool=false,
#														distr_axis::Int=10,
#														kwargs...
#														)::Vector{<:Union{Array,DArray,SharedArray}
#																		} where N 
#
#	@assert 3<=N<=4
#
#	init_overlaps_line(size(wfs,2), nr_kPoints_from_mesh(wfs,dir,N-2);
#								parallel=N==4 && parallel,
#								distr_axis=3-dir,
#								kwargs...
#								)
#
#
#end  
#
#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function init_wlo_mesh(nwf::Int, nmesh::Vararg{Int,N};
#											 kwargs...
#											)::Union{Array{ComplexF64,N+2},
#															 DArray{ComplexF64,N+2},
#															 SharedArray{ComplexF64,N+2},
#															 } where N
#
#	init_storage(ComplexF64, (nwf,nwf), nmesh...; kwargs...)
#
#end 

function init_wlo_mesh_ak(nwf::Int, nmesh::Int...; kwargs...)

	init_storage_ak(ComplexF64, (nwf,nwf), nmesh...; kwargs...)

end 

function init_wlo_mesh_ak(wfs::AbstractArray{ComplexF64,3}, 
											 nmesh::Int=nr_kPoints_from_mesh1(wfs),
											)

	init_wlo_mesh_ak(size(wfs,2), check_nr_kPoints(nmesh))

end   


function init_wlo_mesh_ak(wfs::AbstractArray{ComplexF64,4},
											 nmesh::Int
											 )#::Array{ComplexF64,4}

	init_wlo_mesh_ak(size(wfs,2), check_nr_kPoints(nmesh), nmesh)

end  

function init_wlo_mesh_ak(dir::Int,
													s::NTuple{2,Int},
													n::NTuple{2,Int},
													T...
													;
											 kwargs...
											)

	init_wlo_mesh_ak(s[2], n...; wlo_filter_distr_kwargs(dir; kwargs...)...)

end 

function init_wlo_mesh_ak(dir::Int,
													size_wfmesh::NTuple{4,Int};
											 kwargs...
											)

	(s,n,), = ak_from_storage(size_wfmesh)

	return init_wlo_mesh_ak(dir, s, n; kwargs...)

end


function init_wlo_mesh_ak(dir::Int,
											 wfs::AbstractArray{ComplexF64,4};
											 kwargs...
											)

	init_wlo_mesh_ak(dir, size(wfs); kwargs...) 

end 

function init_wlo_mesh(args...; kwargs...
											 )::Union{Array{ComplexF64},
															 DArray{ComplexF64},
															 SharedArray{ComplexF64}
															 } 

	wrapper_ak(init_wlo_mesh_ak, init_storage, args...; kwargs...)

end 

#function init_wlo_mesh(dir::Int,
#											 wfs::AbstractArray{ComplexF64,4};
#											 parallel::Bool=false,
#											 kwargs...
#											)::Union{Array{ComplexF64,4},
#															 DArray{ComplexF64,4},
#															 SharedArray{ComplexF64,4}
#															 }
#
#	init_wlo_mesh(size(wfs,2), 
#								nr_kPoints_from_mesh(wfs,1),
#								nr_kPoints_from_mesh(wfs,2);
#								parallel=parallel,
#								distr_axis=3-dir,
#								kwargs...
#								)
#end 
#


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function wlo1_on_mesh1_inplace(Bloch_WFs::AbstractArray{ComplexF64,3}
															 )::Array{ComplexF64,3}


	dest = init_wlo_mesh(Bloch_WFs)

	overlaps = init_overlaps_line(Bloch_WFs)

	wlo1_on_mesh1_inplace!(dest, overlaps..., Bloch_WFs)

	return dest 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function wrapper_ak(akfun::Function, 
										funs::Tuple{Vararg{<:Function}},
										args...; kwargs...) 

	a, k = akfun(args...; kwargs...)

	return [f(a...; k...) for f in funs]

end 

function wrapper_ak(akfun::Function, fun::Function, args...; kwargs...) 

	a, k = akfun(args...; kwargs...)

	return fun(a...; k...) 

end 


function wlo1_on_mesh_inplace(dir1::Int, 
															WFs::SubOrSArray{ComplexF64,4};
															kwargs...
															)::SharedArray{ComplexF64,4}

	dest,dest_i = init_storage_ida(dir1, WFs; 
																 custom_ak=init_wlo_mesh_ak,
																 kwargs...,
																 ) 

	overlaps,ov_i = init_overlaps_line_ida(dir1, WFs; kwargs...)

	wlo1_on_mesh_inplace!(dest, overlaps..., WFs, dir1; 
												array_distrib=(dest_i,ov_i))

	return dest 

end 

function wlo1_on_mesh_inplace(dir1::Int, 
															WFs::Union{SubOrArray{ComplexF64,4},
																				 SubOrDArray{ComplexF64,4}
																				 };
															kwargs...
															)::Union{Array{ComplexF64,4},
																			DArray{ComplexF64,4},
																			}

	dest = init_wlo_mesh(dir1, WFs; kwargs...)  

	@assert !isa(dest,SubOrSArray) 

	overlaps = init_overlaps_line(dir1, WFs; kwargs...)


	wlo1_on_mesh_inplace!(dest, overlaps..., WFs, dir1)

	return dest 

end 

function _wlo1_on_mesh_inplace!(
						dest::Td,
						overlaps::To,
						ov_aux::Ta,
						 WFs::Tw,
						 dir1::Int
						)::Nothing where {T<:ComplexF64,
															T4<:Union{SubOrArray{T,4},SubOrSArray{T,4}},
															Td<:T4, Tw<:T4, To<:T4, Ta<:T4}

	@assert size(overlaps)==size(ov_aux) 

	n1 = select_size_dir(overlaps,1)
	n2 = select_size_dir(overlaps,2)
	

	@assert n1==1 || n2==1 "There must be one singleton dimension"

	d = get_mesh_dir(overlaps, n1==n2 ? 3-dir1 : (n1==1 ? 1 : 2))

	return _wlo1_on_mesh_inplace!(dest, 
																selectdim(overlaps,d,1), selectdim(ov_aux,d,1),
																WFs, dir1) 

end  


function _wlo1_on_mesh_inplace!(
						dest::Td,
						overlaps::To,
						ov_aux::Ta,
						 WFs::Tw,
						 dir1::Int
						)::Nothing where {T<:ComplexF64,
															T4<:Union{SubOrArray{T,4},SubOrSArray{T,4}},
															T3<:Union{SubOrArray{T,3},SubOrSArray{T,3}},
															Td<:T4, Tw<:T4, To<:T3, Ta<:T3}


	for d=1:2 
		@assert select_size_dir(WFs,d)==select_size_dir(dest,d)
	end 

	for k2=1:select_size_dir(WFs, 3-dir1) # parameter 

		wlo1_on_mesh1_inplace!(dest, overlaps, ov_aux, WFs, dir1, k2)

	end 

	return 

end 





function wlo1_on_mesh_inplace!(
						 dest::SubOrArray{ComplexF64,4},
						 overlaps::SubOrArray{ComplexF64,3},
						 ov_aux::SubOrArray{ComplexF64,3},
						 WFs::SubOrArray{ComplexF64,4},
						 dir1::Int;
						 array_distrib=nothing,
						)::Nothing

	_wlo1_on_mesh_inplace!(dest, overlaps, ov_aux, WFs, dir1) 

end 


function wlo1_on_mesh_inplace!(
						 dest::SubOrDArray{ComplexF64,4},
						 overlaps::SubOrDArray{ComplexF64,4},
						 ov_aux::SubOrDArray{ComplexF64,4},
						 WFs::AbstractArray{ComplexF64,4}, #### copied on each core!
						 dir1::Int;
						 array_distrib=nothing
						)::Nothing

	map(procs(dest)) do p

		@spawnat p begin      

			_wlo1_on_mesh_inplace!(
											localpart(dest),
											localpart(overlaps),
											localpart(ov_aux),
											localpart_(WFs, dest, 3-dir1),
											dir1, 
											)

		end  

	end .|> fetch 

	return 

end 




function wlo1_on_mesh_inplace!(
						 dest::SubOrSArray{T,4},
						 overlaps::SubOrSArray{T,4},
						 ov_aux::SubOrSArray{T,4},
						 WFs::SubOrSArray{T,4},
						 dir1::Int;
						 array_distrib::Tuple{<:AbstractDict,<:AbstractDict},
						)::Nothing where T<:ComplexF64

	map(procs_(array_distrib)) do p 
	
		li,liov = getindex.(array_distrib,p)

		@spawnat p begin      

			_wlo1_on_mesh_inplace!(
														 localpart_(dest,li),
											localpart_(overlaps,liov), 
											localpart_(ov_aux,liov), 
											localpart_(WFs, li, 3-dir1),
											dir1, 
											)

		end  

	end .|> fetch 

	return 

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function nr_kPoints_from_mesh(W::AbstractArray{<:Number,N})::Int where N

	@assert N>=2 

	@assert size(W,N)==size(W,N-1)

	return size(W,N)+1

end 

function nr_kPoints_from_mesh(W::AbstractVector{<:AbstractArray},
																			args::Int...)::Int 

	only(unique(nr_kPoints_from_mesh(w, args...) for w=W))

end 

function nr_kPoints_from_mesh(W::Union{<:AbstractArray{<:Number},
																			<:Tuple{Vararg{Int}}
																			}, args::Int...)::Int 

	select_size_dir(W, args...) + 1

end 


function nr_kPoints_from_mesh1(W::AbstractArray{<:Number,N})::Int where N

	nr_kPoints_from_mesh(W, 1, 1)

end 

function select_size_dir(A::Tuple{Vararg{Int}}, args...)::Int  

	select_mesh_dir(A, args...)
	
end 

function select_size_dir(A::AbstractArray, args...)::Int  

	select_size_dir(size(A), args...)
	
end 

function select_mesh_dir(I::Tuple{Vararg{T}}, args::Int...
						 )::T where T<:Union{Int, UnitRange{Int}} 

	I[get_mesh_dir(I,args...)]

end 

function get_mesh_dir(N::Int, dir::Int, meshdim::Int=2)::Int 

	@assert N>=meshdim && 1<=dir<=meshdim 

	N-meshdim+dir 

end 

function get_mesh_dir(A::Union{NTuple{N,<:Union{Int,UnitRange{Int}}},
															 AbstractArray{<:Number,N}
															 },
											args::Int...
											)::Int where N 

	get_mesh_dir(N, args...)

end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function procs_(array_distrib::AbstractDict)::Vector{Int}

	k = keys(array_distrib) 

	@assert issubset(k,procs())  

	return sort!(vcat(k...))

end 

function procs_(array_distribs::Tuple{Vararg{<:AbstractDict}}
								)::Vector{Int}

	@assert same_procs(array_distribs)
	
	return procs_(array_distribs[1])

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function localpart_(B::SubOrDArray{<:Number}, 
										)::AbstractArray

	localpart(B) 

end 

function localpart_(B::SubOrDArray{<:Number}, 
										A::SubOrDArray{<:Number}, args... 
										)::AbstractArray

	@assert same_mesh_distrib(B,A) 

	return localpart(B) 

end 

function localpart_(B::Union{<:SubOrArray{<:Number}, <:SubOrSArray{<:Number}},
										I::Tuple{Vararg{UnitRange{Int}}}, 
										dir::Int)::AbstractArray

	d = get_mesh_dir(I,dir) 

	return selectdim(B, d, I[d])

end 


function localpart_(B::SubOrArray{<:Number}, 
										A::SubOrDArray{<:Number}, dir::Int...)::AbstractArray

	localpart_(B, localindices(A), dir...)

end 

function localpart_(B::Union{<:SubOrArray{<:Number}, <:SubOrSArray{<:Number}},
										I::Tuple{Vararg{UnitRange{Int}}}
									 )::AbstractArray

	select_mesh_point(B, select_mesh_dir(I,1), select_mesh_dir(I,2))

end  

function localpart_(F::Function, A::SubOrDArray{<:Number}
#										I::Tuple{Vararg{UnitRange{Int}}}
									 )::Function 

	I1, I2 = get_big_inds(A) 

	function F2(i1::Int, i2::Int) 
		
		F(I1[i1], I2[i2])

	end 

end 
function localpart_(F::Function, 
										A_or_Ainds::Union{<:SubOrDArray{<:Number},
														 <:Tuple{Vararg{UnitRange{Int}}}
														 },
									 )::Function 

	I1, I2 = get_big_inds(A_or_Ainds)

	if A_or_Ainds isa SubOrDArray 
	
		@assert (I1,I2) == get_big_inds(localindices(A_or_Ainds))
	end 


	function F2(i1::Int, i2::Int) 
		
		F(I1[i1], I2[i2])

	end 

end 
																			
function localpart_(A,
										arrays::AbstractVector{<:SubOrDArray},
										args...
										)

	@assert same_mesh_distrib(arrays)

	localpart_(A, arrays[1],args...)

end 

function localpart_(arrays::AbstractVector{<:AbstractArray},
										args...)::Vector{<:AbstractArray}

	[localpart_(a, args...) for a in arrays] 

end 

function localpart_(arrays1::AbstractVector{<:AbstractArray}, 
										arrays2::AbstractVector{<:SubOrDArray}, 
										::Vararg)::Vector{<:AbstractArray}
	
	@assert same_mesh_distrib(arrays2)

	localpart_(arrays1, arrays2[1])

end 


function get_big_inds(
										arrays::AbstractVector{<:SubOrDArray},
										)

	@assert same_mesh_distrib(arrays)
										
	get_big_inds(arrays[1])

end 

function get_big_inds(A::SubOrDArray{<:Number}
											 )

	get_big_inds(localindices(A))

end  

function get_big_inds(I::Tuple{Vararg{UnitRange{Int}}},
										 )

	Tuple(select_mesh_dir(I,d) for d=1:2)

end 


										

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#
#
#function Wannier_subspace_on_mesh_one(ij::NTuple{2,Int},
#																			W::AbstractArray,
#													s::Union{Function,Real}, 
#													)
#
#	Wannier_subspace(get_item_ij(ij, W), s)
#
#end 
#
#
#
#
#function Wannier_subspace_on_mesh(W::AbstractArray,
#													s::Union{Function,Real}, 
#													)
#
#	error() 
#
#	store_on_mesh(Wannier_subspace_on_mesh_one, 
#								nr_kPoints_from_mesh(W), tuple, W, s)
#
#
#end 

#function Wannier_subspaces_on_mesh(W::AbstractArray)
#
#	error() 
#
#	store_on_mesh(Wannier_subspaces∘get_item_ij,
#								nr_kPoints_from_mesh(W), tuple, W) 
#
#end  



function Wannier_eigen_on_mesh!(W::SubOrArray{ComplexF64,4};
																		kwargs... 
																		)::Vector{<:AbstractArray}
	
	@assert !get(kwargs,:parallel,false) 

	_Wannier_eigen_on_mesh!(W; kwargs...)

end 

function Wannier_eigen_on_mesh!(W::Union{SubOrDArray{ComplexF64,4},
																						 SubOrSArray{ComplexF64,4},
																						 };
																		kwargs...
																	 )::Vector{<:AbstractArray}

	@assert get(kwargs,:parallel,false) 

	_Wannier_eigen_on_mesh!(W; wlo_filter_distr_kwargs(;kwargs...)...)

end   


function _Wannier_eigen_on_mesh!(W::Union{SubOrArray{ComplexF64,4},
																							SubOrDArray{ComplexF64,4}
																							};
																		 kwargs...
																		 )::Vector{<:AbstractArray}


	store_on_mesh!!(Wannier_eigen!!, W, 
									init_Wannier_eigen(W; kwargs...))

end 


function _Wannier_eigen_on_mesh!(W::SubOrSArray{ComplexF64,4};
																		 kwargs...
																		)::Vector{SharedArray}


	storage, (storage_i,) = init_storage_ida(init_Wannier_eigen0(W),
																					 nr_kPoints_from_mesh(W);
																					 kwargs...) 

	store_on_mesh!!(Wannier_eigen!!, W, storage; array_distrib=storage_i)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Wannier_eigen_on_mesh!(W::SubOrArray{ComplexF64,3};
#																		 kwargs...
																		 )::Vector{<:AbstractArray}


	store_on_mesh1!!(Wannier_eigen!!, W, init_Wannier_eigen(W))

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#







function Wannier_subspace_on_mesh(s::Union{Function,Real}, Hdata...)

	W,U = wlo1_on_mesh(Hdata...)

	return Wannier_subspace_on_mesh(W, s),W,U

end 


#function Wannier_band_basis_mesh(
#																 W::AbstractArray,
#																 psiH::AbstractArray,
#																 s::Union{Function,Real})
#
#	@warn "Obsolete, output modified"
#
#	eig = Wannier_subspace_on_mesh(W, s)
#
#	return Wannier_band_basis_mesh(eig[1], psiH)[1]
#
#end   

#function Wannier_band_basis_mesh(
#																 eigW1::AbstractArray{<:AbstractArray},
#																 psiH::AbstractArray{ComplexF64,N};
#																 kwargs...
#																 )::AbstractArray{ComplexF64,N} where N
#
#	@assert length(eigW1)==2 
#
#	Wannier_band_basis_mesh(eigW1[1], psiH; kwargs...)
#
#end    
		

function Wannier_band_basis_mesh(
																 Wwf::AbstractArray{ComplexF64,4},
																 psiH::AbstractArray{ComplexF64,4};
																 kwargs...
																 )::AbstractArray{ComplexF64,4}

	store_on_mesh(Wannier_band_basis0, Wannier_band_basis0!,
								nr_kPoints_from_mesh(psiH), tuple, psiH, Wwf; 
								kwargs...)

end    

function Wannier_band_basis_mesh(
																 Wwf::AbstractArray{ComplexF64,3},
																 psiH::AbstractArray{ComplexF64,3};
																 kwargs...
																 )::AbstractArray{ComplexF64,3}

	store_on_mesh1(Wannier_band_basis0, Wannier_band_basis0!,
								nr_kPoints_from_mesh1(psiH), identity, psiH, Wwf; 
								kwargs...)

end  






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function wcc2mesh_fromSubspaces1(dir2::Int,
																 (psiW,nuW)::AbstractVector{<:AbstractArray},
																psiH::AbstractArray{<:Number,4};
																kwargs_...
																)::Vector{<:AbstractArray{Float64,3}}


	kwargs = wlo_filter_distr_kwargs(dir2; kwargs_...)

	wbb_a, = Wannier_band_basis0_ak(psiH, 
																	psi_sorted_energy(psiW; halfspace=true))


	w2,w2_distr = init_storage_ida(dir2, wbb_a...;
																 custom_ak=init_wlo_mesh_ak,
																 kwargs...,
																 ) 

	overlaps,ov_distr = init_overlaps_line_ida(dir2, wbb_a...; kwargs...)

	return wcc2mesh_fromSubspaces1!(w2, overlaps, 
																	dir2, psiW, psiH, 
																	(array_distrib=(w2_distr,ov_distr),); 
																	kwargs...) 

end 

function wcc2mesh_fromSubspaces1!(
					w2, overlaps, 
					dir2, psiW, psiH,
					kwargs_w2;
					kwargs...)

	@assert size(psiW)[1:2] ==(2,2) 


	wbb_distr = inds_distrib_array(psiH, psiW;
																custom_ak=Wannier_band_basis0_ak,
																	 kwargs...) 

#	@warn "psi overwritten, dir2=$dir2"

	store_on_mesh!!(Wannier_band_basis0!, psiW, psiH; array_distrib=wbb_distr) 

	return map([true,false]) do sector  
	
		wbb = psi_sorted_energy(psiH; halfspace=true, occupied=sector) 

		wlo1_on_mesh_inplace!(w2, overlaps..., wbb, dir2; kwargs_w2...)

		return store_on_mesh(get_periodic_eigvals∘Matrix, get_periodic_eigvals!, 
												 w2; kwargs...)

	end  

end 




#
#function wcc2mesh_fromSubspaces1(dir2::Int,
#																 data_dir12::AbstractVector{<:AbstractVector{<:AbstractArray} },
#																psiH::AbstractArray{<:Number,4};
#																kwargs...
#																)::Vector{Array{Float64,3}}
#
#
#	wcc2mesh_fromSubspaces1(dir2, data_dir12[3-dir2], psiH; kwargs...)
#
#end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


	
function wlo2_on_mesh(dir2::Int, data...; 
											distr_axis::Int=3-dir2,
											kwargs...
										 )#::Array{ComplexF64,4}

	wbb = Wannier_band_basis_mesh(data...; distr_axis=3-dir2, kwargs...)

	return wlo1_on_mesh_inplace(dir2, wbb; kwargs...)

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

function get_periodic_eigvals!(dst::AbstractVector{Float64},
															 W::AbstractMatrix)::AbstractVector{Float64}

	copy!(dst, get_periodic_eigvals!(W))

end 


function get_periodic_eigvals!(W::Union{SubOrArray{<:Number,2},SubOrSArray{<:Number,2}},
															 )::Vector{Float64}

	get_periodic_eigvals(LinearAlgebra.eigvals!(W))

end 


function get_periodic_eigvals(W::SubOrArray{<:Number,2}
														 )::Vector{Float64}

	get_periodic_eigvals(LinearAlgebra.eigvals(W))

end  

#function get_periodic_eigvals(W::SubOrDArray{<:Number,2}
#														 )::Vector{Float64}
#
#	get_periodic_eigvals(convert(Array,W))
#
#end 


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

function nuPlusMinus_fromSubspaces((wf_plus,
																		 nu_plus,
																		 wf_minus,
																		 nu_minus,
																		 )::AbstractVector{<:AbstractArray}
									 )::Tuple{<:AbstractArray{<:Real,3},
														<:AbstractArray{<:Real,3}}

	(nu_plus,nu_minus)
	
end 



function polariz_fromSubspaces(storage::AbstractVector{<:AbstractArray}
															 )::Array{Float64,3}

# to be removed 
	sum(nuPlusMinus_fromSubspaces(storage))

end  


function polariz_fromSubspaces1(storage::AbstractVector{<:AbstractArray}
															 )::Array{Float64,3}

	add_nupm(view(storage,[2,4]))


end  

function polariz_fromSubspaces1!(out::AbstractArray{Float64,3},
																 storage::AbstractVector{<:AbstractArray}
															 )::AbstractArray{Float64,3}

	# wf_plus,nu_plus,wf_minus,nu_minus  
	
	add_nupm!(out, view(storage,[2,4]))

end   

function add_nupm(storage::AbstractVector{<:AbstractArray}
															 )::Array{Float64,3}

	sum(storage)

end   

function add_nupm!(out::AbstractArray{Float64,3},
																 storage::AbstractVector{<:AbstractArray}
															 )::AbstractArray{Float64,3}

	copy!(out, storage[1]) # only nu_plus, nu_minus 

	out .+= storage[2] 

	return out 

end   


function polariz_fromSubspaces!(storage::AbstractVector{<:AbstractArray}
									 )::AbstractArray{<:Real,3}


	storage[2] += storage[4] 

	return storage[2] 
	
end 


function WannierGap_fromSubspaces(storage::AbstractVector{<:AbstractArray}
																		)::Float64

	Wannier_min_gap_mesh(nuPlusMinus_fromSubspaces(storage)...)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function check_nu_k_dep(nu::AbstractMatrix{T},
												d::Int;
												atol::Float64=1e-10,
												onewarning::Bool=false,
												assert::Bool=true,
												)::Tuple{Bool,<:AbstractVector{T}} where T<:Real 

	out = true 

	for i=2:size(nu,d) 


		D = Utils.reduce_dist_periodic(max, selectdim(nu,d,i), selectdim(nu,d,1), 1)


		D<atol &&  continue
	
		@warn "nux depends on kx | nuy depends on ky?" 

		@show LinearAlgebra.norm(D)

		@show Utils.reduce_dist_periodic(max, selectdim(nu,3-d,i), selectdim(nu,3-d,1), 1)

		out=false 
#		break 

		onewarning && break 

	end   

	assert && @assert out 

	return (out,selectdim(nu,d,1))

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function Wannier_min_gap_mesh(nu_plus::AbstractArray{<:Real,N},
															nu_minus::AbstractArray{<:Real,N},
															)::Float64 where N

	@assert N>=2 

	out = 1.0

	for i in CartesianIndices(axes(nu_plus)[1:N-2])
		for j in CartesianIndices(axes(nu_minus)[1:N-2])
		
			out = min(out, Utils.reduce_dist_periodic(min, view(nu_plus, i, :, :), view(nu_minus,j, :, :), 1))

		end
	end 

	return out 

end 

function Wannier_min_gap_mesh(nus::AbstractArray{<:Real,N},
															)::Float64 where N


	@assert N>=2  
	
	CI = CartesianIndices(axes(nus)[1:N-2]) 

	out = 1.0 

	for i=1:length(CI)-1

		@views out = min(out, Wannier_min_gap_mesh(nus[CI[i:i],:,:], 
																							 nus[CI[i+1:end],:,:]))

	end 

	return out 

end  



function Wannier_min_gap(nu_plus::Tp,
												 nu_minus::Tm,
												)::Float64 where {
														T<:Union{<:Real,<:AbstractArray{<:Real}},
														Tp<:T, Tm<:T, 
														}

	out = 1.0 

	for i=1:length(nu_plus)

		out = min(out, Utils.reduce_dist_periodic(min, get(nu_plus,i,nu_plus), nu_minus, 1))

	end  

	return out 


end 

function Wannier_min_gap(nu::AbstractVector{<:Real})::Float64 

	out = 1.0 

	for i=2:length(nu)

		out = min(out, Wannier_min_gap(nu[i],view(nu,1:i-1)))

	end 

	return out 

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
function init_Wannier_eigen(W::AbstractArray{<:Number,3};
																kwargs...
																)::Vector{<:AbstractArray} 
 
	init_storage(init_Wannier_eigen0(W), 
							 (nr_kPoints_from_mesh1(W),);
							 kwargs...)

end 

function init_Wannier_eigen(W::AbstractArray{<:Number,4};
																kwargs...
																)::Vector{<:AbstractArray} 
 
	init_storage(init_Wannier_eigen0(W), 
							 nr_kPoints_from_mesh(W);
							 kwargs...)

end 
function init_Wannier_eigen0(W::AbstractArray{<:Number,N},
																)::Vector{<:Array} where N
 

	init_Wannier_eigen0(select_mesh_point(W,Tuple(1 for d=1:N-2)))


end  


function init_Wannier_eigen0(W::AbstractMatrix{t},
																)::Vector{<:Array} where t<:Number 

	n = LinearAlgebra.checksquare(W) 

#	m = div(n,2) # halfspace
	m = n 

	T = promote_type(t,Float64)

	return [
					Matrix{T}(undef, n, m), Vector{Float64}(undef, m),
#					Matrix{T}(undef, n, m), Vector{Float64}(undef, m),
					]

end 


function Wannier_subspaces(W::AbstractMatrix,
													)

	eig = LinearAlgebra.eigen(W)

	E = get_periodic_eigvals(eig.values)

	nu_plus = psien_sorted_energy(eig.vectors, E; 
															halfspace=true, occupied=false)

	nu_minus = psien_sorted_energy(eig.vectors, E; 
															 halfspace=true, occupied=true)

	return vcat(nu_plus,nu_minus)

end 

#function Wannier_eigen!(storage::AbstractVector{<:AbstractArray},
#														W::AbstractMatrix,
#													)
#
#	eig = LinearAlgebra.eigen(W)
#
#	E = get_periodic_eigvals(eig.values)
#
#	eigH!(storage, eig.vectors, E) 
#
#	return 
#
#end 
#


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Wannier_eigen!!(storage::AbstractVector{<:AbstractArray},
														 W::AbstractVector{<:AbstractArray}
														 )::Nothing 

	Wannier_eigen!!(storage, only(W))

end 

function Wannier_eigen!!(storage::AbstractVector{<:AbstractArray},
														 W::Union{<:SubOrSArray{<:Number,2},
																			<:SubOrArray{<:Number,2},
																			},
														 )::Nothing 

	for s in storage 
		@assert isa(s,isa(W,SubOrArray) ? SubOrArray : SubOrSArray)
	end 


	eig = LinearAlgebra.eigen!(W)

	E = get_periodic_eigvals(eig.values)

#	eigH!(view(storage,1:2), eig.vectors, E; halfspace=true, occupied=false)
#	eigH!(view(storage,3:4), eig.vectors, E; halfspace=true, occupied=true)

	eigH!(storage, eig.vectors, E; halfspace=false, rev=true)

	return 

end 

#function sep_Wannier_subspaces(A::AbstractArray{<:Number},
#															 d::Int
#															 )::AbstractVector{<:AbstractArray}
#
#
#	N = size(A,d) 
#
#	u = sortperm_energy(N; halfspace=true,occupied=true) 
#	o = sortperm_energy(N; halfspace=true,occupied=false) 
#
#	return [selectdim(A, 2, o), 
#					selectdim(A, 2, u), 
#					]
#
#
#end 
function sep_Wannier_subspaces((psi,E)::AbstractVector{<:AbstractArray},
															 )::AbstractVector{<:AbstractArray}


	N = size(E,1) 
	@assert size(psi,2)==N 

	o = sortperm_energy(N; halfspace=true,occupied=true) 
	u = sortperm_energy(N; halfspace=true,occupied=false) 

	return [selectdim(psi, 2, o), selectdim(E, 1, o), 

					selectdim(psi, 2, u), selectdim(E, 1, u),
					]


end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#







function Wannier_subspace(W::AbstractMatrix,
													s::Union{Function,Real}, 
													)#::Matrix 

	@warn "Obsolete" 

	eig = LinearAlgebra.eigen(W)

	ov = overlap(eig.vectors)

	occupied_subspace(eig.vectors, 
														 get_periodic_eigvals(eig.values);
														 choose_Wannier_subspace(s)...)



end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_wlo_data_mesh1(psiH::AbstractArray{ComplexF64,3}, 
														occupied::Bool=true
														)::Tuple{Array{ComplexF64,3},
																		#Array{ComplexF64,3},
																		Array{Float64,2},
																		}




	psi = psi_sorted_energy(psiH; halfspace=true, occupied=occupied) 

	W1 = wlo1_on_mesh1_inplace(psi)


	psiW,nuW = Wannier_eigen_on_mesh!(W1)  

	wbb = Wannier_band_basis_mesh(psiW, psi)

	@assert !any(isnan, wbb)

	return (wbb, nuW)

end  



function get_wlo_data_mesh(full_psiH::SubOrArray{ComplexF64,4}, 
									occupied::Bool,
									dir1::Int,
									get_wlo2::Bool;
									kwargs...
									)::Tuple

	psi = psi_sorted_energy(full_psiH; halfspace=true, occupied=occupied) 

	w1 = wlo1_on_mesh_inplace(dir1, psi; kwargs...)


	eigW1 = Wannier_eigen_on_mesh!(w1; dir=dir1, kwargs...)

	do_conv = any(isa(x,SubOrDArray) for x in eigW1) 
	if do_conv 
		eigW1 = convert(Vector{Array}, eigW1)
	end  
	
	eigW1_ret = sep_Wannier_subspaces(eigW1) 

	if get_wlo2 

		wcc2 = wcc2mesh_fromSubspaces1(3-dir1, eigW1, psi; kwargs...)

		return (eigW1_ret, do_conv ? convert(Vector{Array}, wcc2) : wcc2) 

	else 

		return (eigW1_ret, fill(zeros(0,0,0),0))

	end

end 



#
# psi: 0.15-0.18/10_000			 							# j 
# 				even 0.18 
#
# for each dir 
# 	wlo1: 0.12/10_000				 					# k 
# 	subspaces: 0.75/10_000				 			#	j
#		
#		for +/- 
#	 		Wannier basis: 0.5/10_000 				# j,k
#	 		eigvals: 0.015/10_000							# j
# (0.15 + nr_dir * nr_sectors * 2.25) * (n/100)^2 



function get_wlo_data_mesh(full_psiH::SubOrSArray{ComplexF64,4},
									occupied::Bool,
									dir1::Int,
									get_wlo2::Bool;
									kwargs...
									)::Tuple

	psi = psi_sorted_energy(full_psiH; halfspace=true, occupied=occupied) 
	
	W,w_distr1 = init_storage_ida(dir1, psi; 
																 custom_ak=init_wlo_mesh_ak,
																 kwargs...,
																 ) 

	overlaps,ov_distr = init_overlaps_line_ida(dir1, psi; kwargs...)

	wlo1_on_mesh_inplace!(W, overlaps..., psi, dir1; 
												array_distrib=(w_distr1,ov_distr)
												)

	eigW1 = Wannier_eigen_on_mesh!(W; dir=dir1, kwargs...)

	eigW1_ret = sep_Wannier_subspaces(eigW1) 

	get_wlo2 || return (eigW1_ret, fill(zeros(0,0,0),0))


	dir2 = 3-dir1 

	wbb_a, = Wannier_band_basis0_ak(psi, 
																	psi_sorted_energy(eigW1[1]; halfspace=true))

	w_distr2 = inds_distrib_array(dir2, wbb_a...; custom_ak=init_wlo_mesh_ak, kwargs...) 


	wcc2 = wcc2mesh_fromSubspaces1!(Wannier_sector_storage(W),
																	Wannier_sector_storage.(overlaps),
																	dir2, eigW1[1], psi, 
																	(array_distrib=(w_distr2,ov_distr),);
																	kwargs...)
	
	return (eigW1_ret, wcc2) 


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


function Wannier_band_basis0!(
															A::AbstractMatrix{ComplexF64},
															i::Int,
														 psiH::AbstractArray{ComplexF64,3},
														 W1wf_::AbstractArray{ComplexF64,3},
														 )::AbstractMatrix{ComplexF64}

	Wannier_band_basis0!(A, get_item_ij(i, psiH), get_item_ij(i, W1wf_))

end 
 
function Wannier_band_basis0!( 
															psiH::AbstractMatrix{ComplexF64},
														 W1wf_::AbstractMatrix{ComplexF64},
														 )::AbstractMatrix{ComplexF64}

	copy!(psiH, Wannier_band_basis0(psiH, W1wf_))

end 

function Wannier_band_basis0!(
															A::AbstractMatrix{ComplexF64},
															ij::Tuple{Int,Int},
														 psiH::AbstractArray{ComplexF64,4},
														 W1wf_::AbstractArray{ComplexF64,4},
														 )::AbstractMatrix{ComplexF64}

	Wannier_band_basis0!(A, get_item_ij(ij, psiH), get_item_ij(ij, W1wf_))

end 


function Wannier_band_basis0(i::Int,
														 psiH::AbstractArray{ComplexF64,3},
														 W1wf_::AbstractArray{ComplexF64,3},
														 )::Matrix{ComplexF64}

	Wannier_band_basis0(get_item_ij(i, psiH),  get_item_ij(i, W1wf_))

end 

function Wannier_band_basis0(ij::Tuple{Int,Int},
														 psiH::AbstractArray{ComplexF64,4},
														 W1wf_::AbstractArray{ComplexF64,4},
														 )::Matrix{ComplexF64}

	Wannier_band_basis0(get_item_ij(ij, psiH), get_item_ij(ij, W1wf_))

end  

function Wannier_band_basis0(psiH::AbstractMatrix{ComplexF64},
														 W1wf_::AbstractMatrix{ComplexF64},
														 )::Matrix{ComplexF64}

		psiH *  W1wf_

end 

function Wannier_band_basis0!(
															A::AbstractMatrix{ComplexF64},
															psiH::AbstractMatrix{ComplexF64},
														 W1wf_::AbstractMatrix{ComplexF64},
														 )::AbstractMatrix{ComplexF64}

	LinearAlgebra.mul!(A, psiH, W1wf_)

end  

function Wannier_band_basis0_ak(
														 psiH::AbstractArray{ComplexF64,3},
														 W1wf_::AbstractArray{ComplexF64,3},
														 nmesh::Int=nr_kPoints_from_mesh1(psiH);
														 kwargs...
														 ) 

	Wannier_band_basis0_ak(select_mesh_point(psiH,1),
												 select_mesh_point(W1wf_,1),
												 nmesh;
												 kwargs...)

end 

function Wannier_band_basis0_ak(
														 psiH::AbstractArray{ComplexF64,4},
														 W1wf_::AbstractArray{ComplexF64,4},
														 nmesh::Int=nr_kPoints_from_mesh(psiH);
														 kwargs...
														 ) 

	Wannier_band_basis0_ak(select_mesh_point(psiH,1,1),
												 select_mesh_point(W1wf_,1,1),
												 nmesh;
												 kwargs...)

end 

function Wannier_band_basis0_ak(
														 psiH::AbstractArray{ComplexF64,2},
														 W1wf_::AbstractArray{ComplexF64,2},
														 nmesh::Int;
														 kwargs...
														 ) 
	init_storage_ak(ComplexF64, (size(psiH,1),size(W1wf_,2)), (nmesh, nmesh);
								 kwargs...)

end  




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



# w_x in (6.5)

function Wannier_band_basis0(k::AbstractVector, 
														dir1::Int, 
														subspace::Union{Real,Function},
														H::Function,
														Hdata...)
	
	U,E = occupied_subspace(k, H, Hdata...)

#####
#######

	w1 = wlo1(k, dir1, H, Hdata...) 

#######
#
#	W = LinearAlgebra.eigen(w1).vectors 
#
#
##	GramSchmidt!(W,2)
#

#	error()  
########


	W,nu = Wannier_subspace(w1, subspace) 

#####
#
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


	error() 

	psi_H, E, W1, psi_W1, nu, out = Wannier_band_basis0(k, dir1, subspace, H, Hdata...)

####
	@assert nu ≈ get_periodic_eigvals(LinearAlgebra.diag(psi_W1' * W1 * psi_W1))

		eig = LinearAlgebra.eigen(W1)

	@assert only(nu)==get_periodic_eigvals(eig.values[argmax(abs.(overlap(eig.vectors,psi_W1))[:])])

@assert E ≈ LinearAlgebra.diag(psi_H' * H(k, Hdata...) * psi_H)

#	error()

#	println(abs(only(diff(E))))
###



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
