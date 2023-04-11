module CalcWLOadapt
#############################################################################

import Random, LinearAlgebra, Optim


using OrderedCollections: OrderedDict 

import myLibs: ReadWrite, Utils, SignalProcessing

import myLibs.Parameters: UODict 

import ..FILE_STORE_METHOD


import ..MB; MODEL=MB
##import ..BBH; MODEL=BBH 


import ..WLO 
import ..WLO: nr_kPoints, kPoint_start  

import ..Helpers: Symmetries 
#..AdaptiveMesh


import ..CalcWLO
import ..CalcWLO: all_symms_preserved, preserved_symmetries,
									 perturb_strength, 
									 islegend, fn_legend, xxlabel,
									 get_target, FoundFiles, Read,
									 cp_wcc_to_results,
									 get_data_args,
									 set_results!


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

Dependencies = [MODEL,WLO,CalcWLO]

function adaptive_kMesh end 


usedkeys()::Vector{Symbol} = [:kMesh_type, :kMesh_model]

function usedkeys(P::UODict)::Vector{Symbol} 

	adaptive_kMesh(P) ? usedkeys() : Symbol[]

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function adaptive_kMesh(P::UODict)::Bool 

	haskey(P, :kMesh_type) && lowercase(P[:kMesh_type]) == "adaptive"

end 


function kMesh_model(P::UODict)::Function 

	@assert adaptive_kMesh(P)

	n = P[:kMesh_model]

	getf = getproperty(@__MODULE__, Symbol("f_"*n))

	getp = getproperty(@__MODULE__, Symbol(n*"_par_from_vals"))

	return getf∘getp 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function init_results(n::Int, k0::Real,
											obs::AbstractVector{<:AbstractString}
										 )::Dict{String,Any}

	results = CalcWLO.init_results(n, obs) 

# "ks" => WLO.get_kij(n, k0; restricted=false)(1:n)

	@warn "'ks' not yet in 'results'"


	return results

end   







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function symmetrize_and_normalize!(pert::AbstractArray{ComplexF64,4}, 
																	 symms::AbstractString,
																	 args...)::AbstractArray{ComplexF64,4}

	WLO.store_on_mesh!!(Symmetries.symmetrize_HC!, pert)

	if !isempty(MODEL.sepSymmString(symms))

		error("Not ready to symmetrize on non-uniform mesh")

	end 
#	MODEL.symmetrize_on_mesh!(pert, symms, args...)

	WLO.store_on_mesh!!(Symmetries.normalize!, pert)

	return pert 

end 



function get_perturb_on_mesh_(
										 symms::AbstractString, 
										 )::Array{ComplexF64,4}

	symmetrize_and_normalize!(rand(ComplexF64,4,4,n-1,n-1), symms, n, k0) 

end   

function get_perturb_on_mesh_(
										 symms::AbstractString, 
										 n::Int, k0::Real,
										 seed::Int 
										 )::Array{ComplexF64,4}


	Random.seed!(seed) 

	return get_perturb_on_mesh_(symms, n, k0)
 
end 

function get_perturb_on_mesh(symms::AbstractString, n::Int,
														 args...)::Array{ComplexF64,4}
	
	all_symms_preserved(symms) && return zeros(ComplexF64,4,4,n-1,n-1)

	return get_perturb_on_mesh_(symms, n, args...)

end 

function get_perturb_on_mesh(
											P::UODict, args...
										 )::Array{ComplexF64,4}

	get_perturb_on_mesh(preserved_symmetries(P), 
											nr_kPoints(P),
											kPoint_start(P),
											args... 
											) 

end 













#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function get_data(psiH::AbstractArray{ComplexF64,4}, 
									results::AbstractDict;
									parallel::Bool=false
									)::Vector

	error("'get_data' not implemented")

	#(parallel ? pmap : map)(Base.splat(WLO.get_wlo_data_mesh), 
	#												get_data_args(psiH, results))

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function f_line(p::AbstractVector{<:Real}
								)::Function 

	f(x::Real)::Float64 = p[1]*x+p[2] 

end 

function line_par_from_vals(X::AbstractVector{<:Real},
														Y::AbstractVector{<:Real},
														)::Vector{Float64}

#	@assert all(>(0),X)&&all(>(0),Y) "Non-negative values expected"

	@assert !(X[1]≈X[2]) "X values to close from eachother"

	Y[1]≈Y[2] && return zeros(2) 

	@assert !xor(X[1]<X[2],Y[1]<Y[2]) "The function must increase"

	return inv(hcat(X,ones(2)))*Y

end 

function f_square(p::AbstractVector{<:Real}
								)::Function 

	abs2∘f_line(p)

end 

function f_sin(p::AbstractVector{<:Real})::Function 

	l = f_line(view(p,1:2))

	f(x::Real)::Float64 = p[3]*sinpi(l(x))

	return f 

end 

#function f_cos(p::AbstractVector{<:Real})::Function 
#
#	l = f_line(view(p,1:2))
#
#	f(x::Real)::Float64 = p[3]*(1-cospi(l(x)))
#
#	return f 
#
#end 


function f_exp(p::AbstractVector{<:Real})::Function 

	exp∘f_line(p)

end 

function exp_par_from_vals(X::AbstractVector{<:Real},
														Y::AbstractVector{<:Real},
														)::Vector{Float64}

	line_par_from_vals(X,log.(Y))

end 
function sin_par_from_vals(X::AbstractVector{<:Real},
														Y::AbstractVector{<:Real},
														)::Vector{Float64}

	ym,yM = sort(Y)

	return vcat(line_par_from_vals(X, [asin(ym/yM)/pi,0.5]), yM)

end  
#function cos_par_from_vals(X::AbstractVector{<:Real},
#														Y::AbstractVector{<:Real},
#														)::Vector{Float64}
#
#	ym,yM = sort(Y)
#
#	return vcat(line_par_from_vals(X, [1/pi-acos(ym/yM)/pi,0.5]), yM)
#
#end  


function expminv_par_from_vals(X::AbstractVector{<:Real},
														Y::AbstractVector{<:Real},
														)::Vector{Float64}

	line_par_from_vals(X,-inv.(log.(Y)))

end 

function square_par_from_vals(X::AbstractVector{<:Real},
														Y::AbstractVector{<:Real},
														)::Vector{Float64}

	line_par_from_vals(X,sqrt.(Y))


end 

function f_expminv(p::AbstractVector{<:Real})::Function 

	∘(exp, -, inv, f_line(p))

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

bound_rescale_kStep(dk::Real)::Float64 = dk 



function bound_rescale_kStep(dk::Real,
											(dk_min,dk_max)::AbstractVector{<:Real},
											alpha::Real=1
											)::Float64

	max(min(alpha*dk,dk_max),dk_min)

end 


#function smoothmax(a::Real,b::Real)::Float64 
#
#	@assert !xor(a<0, b<0)
#	
#	
#	v = [a,b] 
#
##	sign(argmax(abs, v))
#
#	return sign(a)*LinearAlgebra.norm(v,6)
#
#end 
#
#function smoothmin(a::Real,b::Real)::Float64 
#
#	a+b-smoothmax(a,b) 
#
#end 
#
#softbound(dk::Real, dk_min::Real, dk_max::Real)::Float64 = smoothmin(smoothmax(dk,dk_min),dk_max)
#
#hardbound(dk::Real, dk_min::Real, dk_max::Real)::Float64 = max(min(dk,dk_max),dk_min)
#
#
#
#
#function softBound_rescale_kStep(dk::Real,
#											bounds::AbstractVector{<:Real},
#											alpha::Real=1
#											)::Float64
#
#	softbound(dk*alpha, bounds...)
#
#end 
#
#softBound_rescale_kStep(dk::Real,)::Float64 = dk 
#
#
#function softBound_rescale_kStep(get_dk::Function,
#														 br_args...
#											)::Function 
#
#	function get_dk_(args...)::Float64
#	
#		softBound_rescale_kStep(get_dk(args...), br_args...)
#
#	end 
#
#end 
#softBound_rescale_kStep = bound_rescale_kStep


function bound_rescale_kStep(get_dk::Function,
														 br_args...
											)::Function 

	function get_dk_(args...)::Float64
		
		bound_rescale_kStep(get_dk(args...), br_args...)

	end 

end 
	

function verify_dk_bounds(
											(dk_min,dk_max)::AbstractVector{<:Real},
											N::Int
											)::Vector{Float64}

	K = [dk_min, dk_max]

	if dk_min*(N-1)>2pi 
		@warn "dk_min too large. Using 2pi/(N-1)"

		K[1] = 2pi/(N-1)

	end 

	if dk_max*(N-1)<2pi 

		@warn "dk_max too small. Using 2pi/(N-1)"

		K[2] = 2pi/(N-1)
	end 

	return K

end 

#sum_kSteps(s::Real)::Float64 = s

function sum_kSteps(ks::AbstractVector{<:Real},
										 )::Float64

	abs(ks[end]-ks[1])

end 

function sum_kSteps_dist2pi!(
											 gaps::AbstractVector{Float64},
												 ks::AbstractVector{Float64},
												 args...
													 )::Function 

	function ssi2p(alpha::Union{AbstractVector{<:Real},<:Real})::Float64 
	
		fill_gaps_ks!(gaps, ks, args..., only(alpha))
	
#		println(only(alpha),"\t",sum_kSteps(ks))

		return abs2(2pi-sum_kSteps(ks))

	end 

end  


#
#function rescaling_factor_dk(args...)
#
#	2pi/sum_kSteps(args...)
#
#end 


function init_gaps_ks(nk::Int)::NTuple{2,Vector{Float64}}

	zeros(nk), zeros(nk)

end 



function fill_gaps_ks(nk::Int, 
											args...
											)::NTuple{2,Vector{Float64}}

	fill_gaps_ks!(init_gaps_ks(nk)..., nk, args...)


end 

# merge methods 
function fill_gaps_ks!(
											 gaps::AbstractVector{Float64},
											 ks::AbstractVector{Float64},

												 nk::Int, k0::Real,

												 get_gap_at_k::Union{Function,
																						 SignalProcessing.Dierckx.Spline1D},

											get_dk_for_gap::Function,

											br_args...

											)::Tuple{<:AbstractVector{Float64},
															 <:AbstractVector{Float64}
															 }

	@assert length(gaps)==length(ks)==nk 


	setindex!(ks, 2pi+k0, 1)

	for i=2:nk

		setindex!(gaps, get_gap_at_k(ks[i-1]), i-1)

		setindex!(ks, get_dk_for_gap(gaps[i-1]), i) # ks[i] = temporary dk[i-1]

		setindex!(ks, 
							ks[i-1] - bound_rescale_kStep(ks[i], br_args...),
							i)
# k decreasing  as in WLO 
#

	end 
		
	setindex!(gaps, get_gap_at_k(ks[nk]), nk)

	return gaps, ks

end 



function fill_gaps_ks!(
											 gaps::AbstractVector{Float64},
											 ks::AbstractVector{Float64},

												 nk::Int, k0::Real,

												 (get_gap_at_k,data_gap)::Tuple{Function,<:Any},

											get_dk_for_gap::Function,

											br_args...

											)::Tuple{<:AbstractVector{Float64},
															 <:AbstractVector{Float64}
															 }

	@assert length(gaps)==length(ks)==nk 


	setindex!(ks, 2pi+k0, 1)

	for i=2:nk

		setindex!(gaps, get_gap_at_k(data_gap, ks[i-1]), i-1)

		setindex!(ks, get_dk_for_gap(gaps[i-1]), i) # ks[i] = temporary dk[i-1]

		setindex!(ks, 
							ks[i-1] - bound_rescale_kStep(ks[i], br_args...),
							i)
# k decreasing  as in WLO 
#

	end 
		
	setindex!(gaps, get_gap_at_k(data_gap, ks[nk]), nk)

	return gaps, ks

end 




function find_rescaling_factor_dk(nk::Int, args...) 

	find_rescaling_factor_dk!(init_gaps_ks(nk)..., nk, args...)

end 


function find_rescaling_factor_dk!(
											 gaps::AbstractVector{Float64},
												ks::AbstractVector{Float64}, 

												 nk::Int, k0::Real, 

												 get_gap_at_k,
											get_dk_for_gap::Function,

																	bounds::AbstractVector{<:Real}, 

														 )#::Tuple{Float64,Vector{Float64}}

	@assert length(gaps)==length(ks)==nk  

	bounds_new = verify_dk_bounds(bounds, nk)

	#fill_gaps_ks!(
	#							gaps, ks, 
	#							nk, k0,
	#							get_gap_at_k,
	#							get_dk_for_gap, 
	#							bounds_new,
	#							)

	#@show ks 

	#alpha0 = rescaling_factor_dk(ks)

	alpha0 = 1.0 

#	@show alpha0 

	sksd2p = sum_kSteps_dist2pi!(
								gaps, ks, 
								nk, k0,
								get_gap_at_k,
								get_dk_for_gap, 
								bounds_new,
								)

	sol = Optim.optimize(sksd2p, 
											 #[1e-10], [1e10], 
											 [alpha0])

	alpha = only(Optim.minimizer(sol))

#	@show alpha Optim.minimum(sol)

	get_dk_for_gap_new = bound_rescale_kStep(get_dk_for_gap, bounds_new, alpha)

	return (
					fill_gaps_ks!(gaps, ks, nk, k0, get_gap_at_k, get_dk_for_gap_new),
					get_dk_for_gap_new,
					)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function init_arrays_gap(nk::Int, Hdata...)::Tuple

	WF = WLO.init_storage1(WLO.init_eigH(rand(2), Hdata...)[1], nk)

	WF_occ = WLO.psi_sorted_energy(WF; halfspace=true, occupied=true) 
						# only to init the right size 

	return (WF, WLO.init_overlaps_line(WF_occ))

end 
 

function find_largest_step(ks::AbstractVector{<:Real},
													 )::Tuple{Float64,Int}

	i0::Int = 0 
	dk0::Float64 = 0.0
	dk::Float64 = 0.0

	for i = 2:length(ks)

		dk = abs(ks[i-1]-ks[i]) 

		if dk > dk0

			i0 = i
			dk0 = dk 

		end
	end 

	return dk0,i0 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function calc_gap!(
						 ((WF, overlaps), (occupied,dir1,Hdata), (nk,k0)),
							k2::Float64
							)::Float64

	get_K = WLO.get_kij_constrained(nk, k0, 2, k2, 3-dir1)
	
	WLO.store_on_mesh1!!(WLO.eigH!, nk, get_K, WF, Hdata...)

	WF_occ = WLO.psi_sorted_energy(WF; halfspace=true, occupied=occupied)

	w1 = WLO.wlo1_one_inplace!(overlaps..., WF_occ)

	return WLO.Wannier_min_gap(WLO.get_periodic_eigvals!(w1))

end  

function calc_gap!(data)::NTuple{2,Vector{Float64}}

	nk,k0 = data[3]

	ks = WLO.get_kij(nk,k0)(1:nk-1)

	return [calc_gap!(data, k) for k=ks], ks 

end  





#===========================================================================#
#
# Compute 
#
#---------------------------------------------------------------------------#




function Compute(P::UODict; get_fname=nothing, target=nothing,
								 kwargs...)::Dict  

	Compute_(P, target, get_fname; kwargs...)

end   




function Compute_(P::UODict, target, get_fname::Function;
									kwargs...)::Dict{String,Any}

	obs = get_target(target; kwargs...) 
	
	results = Compute_(P, obs; kwargs...)

	ReadWrite.Write_PhysObs(get_fname(P), FILE_STORE_METHOD, results)

	return isnothing(target) ? results : Utils.dict_keepkeys(results, obs) 

end 





function Compute_(P::UODict, target, get_fname::Nothing=nothing; 
										kwargs...)::Dict{String,Any}
	
	
	@assert adaptive_kMesh(P) # otherwise CalcWLO.Compute 



	nk = nr_kPoints(P)
	k0 = kPoint_start(P)
	symms = preserved_symmetries(P)
	strength = perturb_strength(P)


	results = init_results(nk, k0, get_target(target; kwargs...))

	model = kMesh_model(P)


#	results = Dict{String,Any}()
	
#	@show results
#
#
#	perturb1 = get_perturb_on_mesh(P)
#
#	@show LinearAlgebra.norm(perturb1)



#	nk_unif = max(5,div(nk,3))

#	psi = MODEL.get_psiH(P, nk, k0, perturb1, strength)


#	set_results!(results, nk, k0, get_data(psi, results))



## Hamiltonian from P 

	Hdata = MODEL.get_pertHdata(MODEL.params_fromP(P), MODEL.H);

###  
# choose a sector 

	dir1 = 1 
	occupied = true 

		extr_gaps = [0.02, 0.5]
		extr_dk = [1e-7, pi/10]

#### initialize storage  



#	arrays_gap = init_arrays_gap(nk, Hdata...) 



#	for dir1=1:2 

		### pack data  
		
		Hdata_gap = (occupied, dir1, Hdata)

#		data = (arrays_gap, Hdata_gap, (nk,k0))  #
	
nks = [17,71,1000]

		data_1 = (init_arrays_gap(nks[1], Hdata...),
							Hdata_gap, 
							(nks[1], k0)
							)

		gaps_1, ks_1 = calc_gap!(data_1)
		
		extr_gaps[1] = minimum(gaps_1) 

		@show extr_gaps 

		@show sum_kSteps(ks_1)


		dk_from_gap_1 = model(extr_gaps, extr_dk)

		data_2 = (init_arrays_gap(nks[2], Hdata...), Hdata_gap, (nks[2],k0)) 

		(gaps_2, ks_2),dk_from_gap_2 = find_rescaling_factor_dk(nks[2], k0, 
																														(calc_gap!,data_2),
																								dk_from_gap_1, extr_dk)


		reverse!(ks_2)
		reverse!(gaps_2)


		interp_gap = SignalProcessing.Interp1D(ks_2, gaps_2, 3) 
		@show interp_gap isa Function 

		@show typeof(interp_gap)


		(gaps_3, ks_3),dk_from_gap_3 = find_rescaling_factor_dk(nks[3], k0, 
																														interp_gap,
																								dk_from_gap_2, extr_dk)


		dk = sum_kSteps(ks_3)

		eps = ks_3[1]>ks_3[end] ? dk-2pi : 2pi-dk

		dkmax,imax = find_largest_step(ks_3)

		dkmax > 10eps || @warn string("Large deviation from 2pi ",dkmax/eps)


		ks_3[imax:end] .+= eps 

		@show extrema(abs.(diff(ks_3)))

		@assert sum_kSteps(ks_3)≈2pi 
		


#	end 

### 


results["data"]=[nk,WLO.get_kij(nk,k0;restricted=false),
								 (init_arrays_gap(nk, Hdata...), Hdata_gap, (nk,k0)),
								 data_2,
								 dk_from_gap_1
								 ]









	return results

end 







































































































#############################################################################
end # module CalcWLO
