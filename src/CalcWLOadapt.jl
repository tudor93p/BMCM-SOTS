module CalcWLOadapt
#############################################################################

#import Random, LinearAlgebra, Optim
#
#
#using OrderedCollections: OrderedDict 
#
#import myLibs: ReadWrite, Utils, SignalProcessing
#
#import myLibs.Parameters: UODict 
#
#import ..FILE_STORE_METHOD
#
#
#import ..MB; MODEL=MB
###import ..BBH; MODEL=BBH 
#
#
#import ..WLO 
#import ..WLO: nr_kPoints, kPoint_start  
#
#import ..Helpers: Symmetries 
##..AdaptiveMesh
#
#
#import ..CalcWLO
#import ..CalcWLO: xxlabel, get_target, FoundFiles, Read
#
#
##===========================================================================#
##
##
##
##---------------------------------------------------------------------------#
#
#Dependencies = [MODEL,WLO,CalcWLO]
#
#function adaptive_kMesh end 
#
#
#usedkeys()::Vector{Symbol} = [:kMesh_type, :kMesh_model]
#
#function usedkeys(P::UODict)::Vector{Symbol} 
#
#	adaptive_kMesh(P) ? usedkeys() : Symbol[]
#
#end 
#
##===========================================================================#
##
##
##
##---------------------------------------------------------------------------#
#
#
#
#
#function adaptive_kMesh(P::UODict)::Bool 
#
#	haskey(P, :kMesh_type) && lowercase(P[:kMesh_type]) == "adaptive"
#
#end 
#
#
#function kMesh_model(P::UODict)::Function 
#
#	@assert adaptive_kMesh(P)
#
#	n = P[:kMesh_model]
#
#	getf = getproperty(@__MODULE__, Symbol("f_"*n))
#
#	getp = getproperty(@__MODULE__, Symbol(n*"_par_from_vals"))
#
#	return getf∘getp 
#
#end 
#


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function init_results(n::Int,
#											(kx,ky)::AbstractVector{<:AbstractVector{<:Real}},
#											obs::AbstractVector{<:AbstractString}
#										 )::Dict{String,Any}
#
#	results = CalcWLO.init_results(n, obs) 
#
#	results["ks"] = hcat(ky,kx) # consistent: klabels=["k_y","k_x"],W=zeros(n,2)
#
#	@assert issubset(xxlabel(), keys(results))
#
#	return results
#
#end   







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




##===========================================================================#
##
##
##
##---------------------------------------------------------------------------#
#
#function f_line(p::AbstractVector{<:Real}
#								)::Function 
#
#	f(x::Real)::Float64 = p[1]*x+p[2] 
#
#end 
#
#function line_par_from_vals(X::AbstractVector{<:Real},
#														Y::AbstractVector{<:Real},
#														)::Vector{Float64}
#
##	@assert all(>(0),X)&&all(>(0),Y) "Non-negative values expected"
#
#	@assert !(X[1]≈X[2]) "X values to close from eachother"
#
#	Y[1]≈Y[2] && return zeros(2) 
#
#	@assert !xor(X[1]<X[2],Y[1]<Y[2]) "The function must increase"
#
#	return inv(hcat(X,ones(2)))*Y
#
#end 
#
#function f_square(p::AbstractVector{<:Real}
#								)::Function 
#
#	abs2∘f_line(p)
#
#end 
#
#function f_sin(p::AbstractVector{<:Real})::Function 
#
#	l = f_line(view(p,1:2))
#
#	f(x::Real)::Float64 = p[3]*sinpi(l(x))
#
#	return f 
#
#end 
#
##function f_cos(p::AbstractVector{<:Real})::Function 
##
##	l = f_line(view(p,1:2))
##
##	f(x::Real)::Float64 = p[3]*(1-cospi(l(x)))
##
##	return f 
##
##end 
#
#
#function f_exp(p::AbstractVector{<:Real})::Function 
#
#	exp∘f_line(p)
#
#end 
#
#function exp_par_from_vals(X::AbstractVector{<:Real},
#														Y::AbstractVector{<:Real},
#														)::Vector{Float64}
#
#	line_par_from_vals(X,log.(Y))
#
#end 
#function sin_par_from_vals(X::AbstractVector{<:Real},
#														Y::AbstractVector{<:Real},
#														)::Vector{Float64}
#
#	ym,yM = sort(Y)
#
#	return vcat(line_par_from_vals(X, [asin(ym/yM)/pi,0.5]), yM)
#
#end  
##function cos_par_from_vals(X::AbstractVector{<:Real},
##														Y::AbstractVector{<:Real},
##														)::Vector{Float64}
##
##	ym,yM = sort(Y)
##
##	return vcat(line_par_from_vals(X, [1/pi-acos(ym/yM)/pi,0.5]), yM)
##
##end  
#
#
#function expminv_par_from_vals(X::AbstractVector{<:Real},
#														Y::AbstractVector{<:Real},
#														)::Vector{Float64}
#
#	line_par_from_vals(X,-inv.(log.(Y)))
#
#end 
#
#function square_par_from_vals(X::AbstractVector{<:Real},
#														Y::AbstractVector{<:Real},
#														)::Vector{Float64}
#
#	line_par_from_vals(X,sqrt.(Y))
#
#
#end 
#
#function f_expminv(p::AbstractVector{<:Real})::Function 
#
#	∘(exp, -, inv, f_line(p))
#
#end 
#
#
#
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



#sum_kSteps(s::Real)::Float64 = s
#
#function sum_kSteps(ks::AbstractVector{<:Real},
#										nk::Int
#										 )::Float64
#
#	check_vec_len(nk,ks)
#
#	@assert length(ks)==nk 
#
#	return abs(ks[end]-ks[1])  
#
#end 
#
#function sum_kSteps_dist2pi!(
#											 gaps::AbstractVector{Float64},
#											 ks::AbstractVector{Float64},
#												 nk::Int, 
#												 args...
#													 )::Function 
#
#	function ssi2p(alpha::Union{AbstractVector{<:Real},<:Real})::Float64 
#	
#		fill_gaps_ks!(gaps, ks, nk, args..., only(alpha))
#	
##		println(only(alpha),"\t",sum_kSteps(ks))
#
#		return abs(2pi-sum_kSteps(ks,nk))
#
#	end 
#
#end  
#
#
#
#
#function init_gaps_ks(nk::Int)::NTuple{2,Vector{Float64}}
#
#	Tuple(WLO.init_storage(Float64, nk+1) for i=1:2)
#
#end 
#
#function check_vec_len(nk::Int, vectors::AbstractVector{<:Real}...)
#
#	for v in vectors 
#		@assert WLO.nr_kPoints_from_mesh1(v)==nk+1
#	end 
#
#end 
#
#
#
#
#function fill_gaps_ks(nk::Int, 
#											args...
#											)::NTuple{2,Vector{Float64}}
#
#	fill_gaps_ks!(init_gaps_ks(nk)..., nk, args...)
#
#
#end 
#
#
#function next_kPoint(get_dk_for_gap::Function,
#										 gap::Real,
#										 k_prev::Real,
#										br_args...)
#
#	dk = get_dk_for_gap(gap)
#
#	@assert dk>0 "Decrease minimum expected gap"
#
#	WLO.next_kPoint(k_prev, 
#									bound_rescale_kStep(dk, br_args...)
#									)
#
#end 
#
#
##===========================================================================#
##
##
##
##---------------------------------------------------------------------------#
#
#function wrapper_get_gap(get_gap_at_k::Union{Function,
#																						 SignalProcessing.Dierckx.Spline1D},
#												 k::Real
#												 )::Float64
#
#
#	get_gap_at_k(k)
#
#end 
#function wrapper_get_gap(
#												 (get_gap_at_k,data_gap)::Tuple{Function,<:Any},
#												 k::Real
#												 )::Float64
#
#	get_gap_at_k(data_gap, k)
#
#end 
#
#
#function fill_gaps_ks!(
#											 gaps_::AbstractVector{Float64},
#											 ks_::AbstractVector{Float64},
#
#												 nk::Int, k0::Real,
#																	unif_kij::Function,
#																	uniq_kinds::AbstractVector{Int},
#																	ind_minusk::Function,
#																	
#																	get_gap_at_k::Union{<:Tuple{Function,<:Tuple}, Function, SignalProcessing.Dierckx.Spline1D},
#
#											get_dk_for_gap::Function,
#
#											br_args...
#
#											)::Tuple{<:AbstractVector{Float64},
#															 <:AbstractVector{Float64}
#															 }
#
#	check_vec_len(nk, gaps_, ks_)
#
#	gaps = view(gaps_, uniq_kinds)
#	ks = view(ks_, uniq_kinds)
#
#	setindex!(ks, unif_kij(uniq_kinds[1])[1], 1)
#
#	for i=2:length(uniq_kinds)
#
#		setindex!(gaps, wrapper_get_gap(get_gap_at_k, ks[i-1]), i-1)
#
#		setindex!(ks, 
#							next_kPoint(get_dk_for_gap, gaps[i-1], ks[i-1], br_args...),
#							i) 
#
#	end  
#		
#	setindex!(gaps, wrapper_get_gap(get_gap_at_k, ks[end]), length(gaps)) 
#
#	for i in uniq_kinds 
#
#		j = ind_minusk(i)
#
#		setindex!(ks_, -ks_[i], j)
#
#		setindex!(gaps_, gaps_[i], j)
#
#	end 
#	
#	
#	#	ks_[nk] and gaps_[nk] unaccessed so far 
#		
#	setindex!(ks_, 
#						next_kPoint(get_dk_for_gap, gaps_[nk-1], ks_[nk-1], br_args...), 
#						nk)
#
#	setindex!(gaps_, wrapper_get_gap(get_gap_at_k, ks_[nk]), nk)
#
#
#
#	return gaps, ks
#
#end 
#
#
#
#
#function find_rescaling_factor_dk(nk::Int, args...; kwargs...)
#
#	find_rescaling_factor_dk!(init_gaps_ks(nk)..., nk, args...; kwargs...)
#
#end 
#
#
#function find_rescaling_factor_dk!(
#											 gaps::AbstractVector{Float64},
#												ks::AbstractVector{Float64}, 
#
#												 nk::Int, k0::Real, 
#
#												 get_gap_at_k,
#
#											get_dk_for_gap::Function,
#
#																	bounds::AbstractVector{<:Real};
#																	optim_tol::Float64=1e-8
#
#								
#
#														 )#::Tuple{Float64,Vector{Float64}}
#
##@assert length(gaps)==length(ks)==nk  
#	check_vec_len(nk, gaps, ks)
#
#	bounds_new = verify_dk_bounds(bounds, nk) 
#
#	uniq_kinds = WLO.uniqueInds_kMirror(nk,k0)
#	unif_kij = WLO.get_kij(nk,k0) #restricted=false)
#	ind_minusk = WLO.ind_minusk(nk,k0)
#
#
#	fill_gaps_ks!(gaps, ks, nk, k0, 
#																	unif_kij,
#																	uniq_kinds,
#																	ind_minusk,
#
#								get_gap_at_k, get_dk_for_gap, 
#								bounds_new,
#																	)
#
#
#	sksd2p = sum_kSteps_dist2pi!(
#								gaps, ks, 
#								nk, k0,
#																	unif_kij,
#																	uniq_kinds,
#																	ind_minusk,
#
#
#								get_gap_at_k,
#								get_dk_for_gap, 
#								bounds_new,
#								)
#	sol = Optim.optimize(sksd2p, [1.0],
#											 Optim.Options(g_tol = optim_tol,
#																		 f_abstol=optim_tol,
#																		 ))
#
#	alpha = only(Optim.minimizer(sol))
#
#	get_dk_for_gap_new = bound_rescale_kStep(get_dk_for_gap, bounds_new, alpha)
#
#	fill_gaps_ks!(gaps, ks, nk, k0, 
#																	unif_kij,
#																	uniq_kinds,
#																	ind_minusk,
#								get_gap_at_k, get_dk_for_gap_new)
#
#
#
##	correct_sum_kSteps!(ks, nk, uniq_kinds, ind_minusk)
#
#	return (gaps, ks), get_dk_for_gap_new 
#
#end 
#
#

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




#function correct_sum_kSteps!(ks::AbstractVector{Float64},
#														 nk::Int,
#														 uniq_kinds::AbstractVector{Int},
#														 ind_minusk::Function,
#														 )::AbstractVector{Float64} 
#
#	dk = sum_kSteps(ks,nk)
#	
#	#	ks[end] == next_kPoint(ks[1],2pi) = ks[1]-2pi 
#
#	#eps = ks[1]>ks[end] ? dk-2pi : 2pi-dk
#
#	eps = abs(dk-2pi)
#
#	dkmax,imax = find_largest_step(view(ks,uniq_kinds))
#
#
#	@assert dkmax > 10eps string("Large deviation from 2pi: dk=$dkmax, eps=$eps")
#
#	i = uniq_kinds[imax]
#
#	dk = diff(ks[max(1,i-1):min(i+1,end)])[1] 
#
#	i_zero = argmin(Utils.dist_periodic(ks,0,2pi))
#	i_pi = argmin(Utils.dist_periodic(ks,pi,2pi))
#
#
#	dk>0 # increasing 
#	dk<0 # decreasing 
#
## i0 invariant to mirror might exist 
## 
#
#
#error() 
#
#	return ks 
#
#end 
#
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


#function calc_gap!(data)::NTuple{2,Vector{Float64}}
#
#	nk,k0 = data[2]
#
#	ks = WLO.get_kij(nk,k0)(1:nk-1)
#
#	return [calc_gap!(data, k) for k=ks], ks 
#
#end  
#



##===========================================================================#
##
##
##
##---------------------------------------------------------------------------#
#
#"""
#sector = (occupied,dir1)
#nks = (rough_nk, decent_nk, dense_nk)
#
#Compute the gap of W_dir1 for k_dir1 uniform and k_dir2 adaptive.
#
#Find the parameter vector ks_dir2 which is more dense for small gaps
#"""
#function threestep_find_ks(Hdata, (nks,k0), 
#													 sector,
#						(extrema_gaps, extrema_dk),
#						model::Function
#						)
#
#	@assert length(nks)==3 && issorted(nks)
#
#
#	# --- rough estimation of gap(k) in a halfspace --- #
#	
#	data_1 = pack_data_gap(Hdata, (nks[1],k0), sector)
#
#	gaps_1, ks_1 = [calc_gap!(data_1,k) for k=WLO.uniqueKs_kMirror(nks[1],k0)]
#
#
#
#	min_gap = min(extrema_gaps[1], minimum(gaps_1))
#	
#	dk_from_gap_1 = model([min_gap, extrema_gaps[2]], extrema_dk)
#	
#
#
#
#	# --- linearly scale 'model' such that ks span a halfspace --- #
#
#	gap_at_k_2 = (calc_gap!,pack_data_gap(Hdata, (nks[2],k0), sector)) 
#
#	(gaps_2, ks_2),dk_from_gap_2 = find_rescaling_factor_dk(
#																													nks[2], k0, 
#																													gap_at_k_2,
#																							dk_from_gap_1, extrema_dk;
#																							optim_tol=1e-8
#																							)
#
#
#	# --- scale again, full precision --- # 
#
#	if ks_2[2]<ks_2[1] 
#
#		reverse!(ks_2)
#		reverse!(gaps_2)
#
#	end 
#
#	gap_at_k_3 = SignalProcessing.Interp1D(ks_2, gaps_2, 3) 
#	
#	(gaps_3, ks_3),dk_from_gap_3 = find_rescaling_factor_dk(nks[3], k0, 
#																													gap_at_k_3,
#																							dk_from_gap_2, extrema_dk;
#																							optim_tol=1e-14
#																						)
#
#	return (gaps_3, ks_3, dk_from_gap_3,
#					(nks[3],k0),
#					sector)
#
#end 	
#
#
#function kxy_from_threestep(
#														 (gaps_1, precomp_ks_1, get_dk_1,
#															(nk_1,k0_1), (occ_1, perp_dir_1)
#															),
#														 (gaps_2, precomp_ks_2, get_dk_2,
#															(nk_2,k0_2), (occ_2, perp_dir_2)
#															),
#														 ;
#														 kwargs...
#														 )::Vector{Vector{Float64}}
#
#
#	@assert perp_dir_1!=perp_dir_2 
#
#	perp_dir_1==2 ? [precomp_ks_1,precomp_ks_2] : [precomp_ks_2,precomp_ks_1]
#
#end 
#





#function Compute_(P::UODict, target, get_fname::Nothing=nothing; 
#										kwargs...)::Dict{String,Any}
#	
#	
#	adaptive_kMesh(P) || return CalcWLO.Compute_(P, target, get_fname;
#																							 kwargs...)
#
#
#	nk = nr_kPoints(P) # same 
#	kxy = calc_kxy_adaptive(P) # diff 
#	results = init_results(nk, kxy, get_target(target; kwargs...))  #diff 
#
#
#
#	return CalcWLO.Compute_!(results, P, kxy) 
#
#
#end 
#






































































































#############################################################################
end # module CalcWLO

