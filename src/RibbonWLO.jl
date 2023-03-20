module RibbonWLO  
#############################################################################

import Statistics, LinearAlgebra, Random

using OrderedCollections: OrderedDict 

import myLibs: ReadWrite, Utils, Lattices, TBmodel

import myLibs.Parameters: UODict 

import ..FILE_STORE_METHOD

import ..WLO, ..Helpers#, ..Symmetries

import ..CalcWLO 
import ..CalcWLO: MODEL 

import ..WLO: nr_kPoints, kPoint_start 

import ..CalcWLO: all_symms_preserved, fn_legend, fn_nolegend, add_legend!

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

Dependencies = [CalcWLO]


usedkeys::Vector{Symbol} = [:width]
									 

calc_observables::Vector{String} = ["WannierBands1",
																		#"Polarization",
																		"WannierDensity"]

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function width(P::UODict)::Int 

	P[:width]

end 



function xxlabel()::Vector{String}

	["x","xlabel"]

end   

#
#isauxfile = isxxlabel = in(xxlabel()) 

function xxlabel(data::AbstractDict)::Tuple{Vector{Int}, String}

	Tuple(getindex(data,k) for k=xxlabel())

end 




prep_obs(::Nothing=nothing)::Vector{String} = String[]

function prep_obs(obs::AbstractString)::Vector{String} 

	obs1 = fn_nolegend(obs)

	obs1 in calc_observables || return prep_obs()

	return vcat(obs1, fn_legend(obs1), xxlabel()) 

end  


function prep_obs(obs::AbstractVector{<:AbstractString})::Vector{String} 

	mapreduce(prep_obs, vcat, obs; init=prep_obs())

end 

function prep_obs(args...)::Vector{String} 
	
	mapreduce(prep_obs, vcat, args; init=prep_obs())

end 



get_target = prep_obs ∘ Helpers.f_get_target(:observables)	   


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

function init_results(w::Int, 
											obs::AbstractVector{<:AbstractString}
										 )::Dict{String,Any}


	results = Dict{String,Any}("x" => Vector(1:2w),
														 "xlabel" => "Eigenvalue number",
														 )

	@assert Set(keys(results))==Set(xxlabel())

	for a in ("WannierBands1",)

		in(a,obs) || continue 

		results[a] = zeros(2w,2)

		add_legend!(results, a, ["dir"])

	end  
	
#	for a in ("Polarization",)
#
#		in(a,obs) || continue 
#
#		results[a] = zeros(w,2)
#
#		add_legend!(results, a, ["dir"])
#
#	end  

	for a in ("WannierDensity",)

		in(a,obs) || continue 

		results[a] = zeros(w,2w,2)

		add_legend!(results, a, ["dir"])

	end  



	return results

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


#function get_data_args(psiH::AbstractArray{<:Number,4},
#											 )::Vector{Tuple{Array,Bool,Int,Bool}}
#
#	[(psiH, true, 1, false), (psiH, true, 2, false), ]
#
#end  
#
#function get_data(psiH::AbstractArray{ComplexF64,4};
#									parallel::Bool=false
#									)::Vector
#
#	(parallel ? pmap : map)(Base.splat(WLO.get_wlo_data_mesh), 
#													get_data_args(psiH))
#
#end 




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

		D = Utils.reduce_dist_periodic(max, selectdim(nu,d,i), circshift(selectdim(nu,d,1),1), 1)

		D<atol &&  continue  

		D = Utils.reduce_dist_periodic(max, selectdim(nu,d,i), circshift(selectdim(nu,d,1),-1),1)

		D<atol &&  continue  

		@warn "nux depends on kx | nuy depends on ky?" 

		@show LinearAlgebra.norm(D)

		out=false 

		onewarning && break 

	end   

	assert && @assert out 

	return (out,selectdim(nu,d,1))

end  


function cp_wcc_to_results(results::AbstractDict{<:AbstractString,<:Any},
													obs::AbstractString,
													wccs::AbstractMatrix{<:Real},
													dir1::Int,
													) 

	copy!(selectdim(results[obs],2,dir1), 
				check_nu_k_dep(wccs, 2)[2]
				)

end  


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function set_results!(results::AbstractDict,
#													nk::Int, k0::Real,
#													data)::Nothing
#
#	for (dir1,d) in enumerate(data)
#					
#		set_results_onedir!(results, nk, k0, dir1, d) 
#
#	end 
#
#end  
#
#function hWpd(
#							 i_atom::Int, 
#							 nr_orb::Int,
#								wbb::AbstractMatrix{ComplexF64},
#											 )
#
#	I = TBmodel.Hamilt_indices(1:nr_orb, i_atom, nr_orb)  
#
#	return sum(Statistics.mean(abs2, selectdim(wbb,1,I), dims=2))
#
#end 

function Wannier_density(
											 wbb::AbstractArray{ComplexF64,3},
											 nr_at::Int
											 )::Matrix{Float64}

	Wannier_density!(zeros(nr_at, size(wbb,2)), wbb) 

end 

function Wannier_density!(
											dest::AbstractMatrix{Float64},
											 wbb::AbstractArray{ComplexF64,3},
											 )::AbstractMatrix{Float64}

	nr_at = size(dest,1)

	nr_k = size(wbb,3)+1 

	nr_wf = size(wbb,2) 

	@assert size(dest,2)==nr_wf 

	nr_orb, aux = divrem(size(wbb,1), nr_at) 

	@assert aux==0 

	rho0 = Statistics.mean(abs2, wbb, dims=3) 

	dest .= 0  

	for orbital=1:nr_orb 

		dest .+= view(rho0, TBmodel.Hamilt_indices(orbital, 1:nr_at, nr_orb),:,1)

	end 

	return dest 

end 

function polarization(rho::AbstractMatrix{<:Real},
											 nu_val::AbstractVector{<:Real},
											 )::AbstractVector{Float64}


#	nu = Utils.bring_periodic_to_interval.(nu_val, 0, 0.99999, 1)

#	nu =nu_val

	rho*nu_val 

end 

function polarization!(dest::AbstractVector{Float64},
											 rho::AbstractMatrix{<:Real},
											 nu_val::AbstractVector{<:Real},
											 )::AbstractVector{Float64}


#	nu = Utils.bring_periodic_to_interval.(nu_val, 0, 0.99999, 1)

#	nu =nu_val

	LinearAlgebra.mul!(dest, rho, nu_val) 

	for f in (isnan,isinf,ismissing), p in dest 

			if f(p) 

				@show dest findall(f,wbb) findall(f,rho0) nu_val 

				error() 

			break  

			end 

	end 

	return dest 

end 


function polarization!(dest::AbstractVector{Float64},
											 wbb::AbstractArray{ComplexF64,3},
											 nu_val::AbstractVector{<:Real}
											 )

	@assert length(nu_val)==size(wbb,2) 

	polarization!(dest, Wannier_density(wbb, length(dest)), nu_val)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function get_nr_at(results::AbstractDict)::Int 

	@assert any(in(calc_observables),keys(results))

#	haskey(results, "Polarization") && return size(results["Polarization"],1) 

	haskey(results, "WannierDensity") && return size(results["WannierDensity"],1)

	haskey(results, "WannierBands1") && return div(size(results["WannierBands1"],1),2)


end 


function set_results_onedir!(results::AbstractDict, 
													dir1::Int,
												 (wbb,nu_val)::Tuple{
														#				<:AbstractArray{ComplexF64,3},
																		<:AbstractArray{ComplexF64,3},
																		<:AbstractArray{Float64,2}}
													)::Nothing

	nu = check_nu_k_dep(nu_val, 2)[2]

	rho = Wannier_density(wbb, get_nr_at(results))


	if haskey(results, "WannierBands1") 

		copy!(selectdim(results["WannierBands1"],2,dir1), nu)

#		cp_wcc_to_results(results, "WannierBands1", nu_val, dir1)

	end 

	if haskey(results, "WannierDensity")

		copy!(selectdim(results["WannierDensity"],3,dir1), rho)

	end 


#	if haskey(results, "Polarization")
#
#		polarization!(selectdim(results["Polarization"],2,dir1), rho, nu)
#
#	end 


	return 

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


function lattice(w::Int, dir2::Int)::Lattices.Lattice 

	latt = Lattices.SquareLattice()

	Lattices.Superlattice!(latt, setindex!(ones(Int,2), w, dir2))

	Lattices.ReduceDim!(latt, dir2)
	
end 

function add_perturb_wrapper(h0::Function)::Function 

	function h(k::AbstractVector{<:Real})::Matrix{ComplexF64}
		
		h0(k)

	end 

	function h(i::Int, get_k::Function, perturb::AbstractArray{<:Number, 3}
						 )::Matrix{ComplexF64}


		out = h(WLO.get_item_ij(i, get_k)) 

#		println(LinearAlgebra.norm(out), " ",LinearAlgebra.norm(WLO.get_item_ij(i, perturb)))

		out .+= WLO.get_item_ij(i, perturb)

		return out 

	end 

end 


function Bloch_Hamilt(w::Int, dir1::Int,#latt::Lattices.Lattice,
											P::UODict)::Function 
	
	latt = lattice(w, 3-dir1)

	l3 = Lattices.NearbyUCs(latt)

	h0 = TBmodel.Bloch_Hamilt(l3; nr_orb=4, Hopping=MODEL.get_hoppf(P))

	return add_perturb_wrapper(h0) 

end 

function Compute_(P::UODict, target, get_fname::Nothing=nothing; 
										kwargs...)::Dict{String,Any}

	nk = nr_kPoints(P) 

	w = width(P) 

	k0 = kPoint_start(P)

#	results = init_results(w, get_target(target; kwargs...))
	results = init_results(w, get_target(calc_observables; kwargs...))


	@assert all_symms_preserved(P)




############### test  

#	perturb2 = CalcWLO.get_perturb_on_mesh("Ct", nk, k0, 1993)  

	Random.seed!(1993)

	s  = 1e-10 

	perturb0 = rand(ComplexF64,4,4) 
	perturb0 .+= perturb0' 

	perturb1 = rand(ComplexF64,4,4) 

	LinearAlgebra.normalize!(perturb0)
	LinearAlgebra.normalize!(perturb1)

	perturb0 .*= s 
	perturb1 .*= s 

	perturb = zeros(ComplexF64, 4w, 4w, nk-1)

	for j=1:w 

		J = TBmodel.Hamilt_indices(1:4,j,4)

		for k=axes(perturb,3) 

			setindex!(perturb, perturb0, J, J, k)

			if j<w 
				setindex!(perturb, perturb1, 
								J, TBmodel.Hamilt_indices(1:4,j+1,4), k)

				end 
			if j>1

				setindex!(perturb, perturb1', 
								J, TBmodel.Hamilt_indices(1:4,j-1,4), k)

			end 

		end 

	end 

#	@show perturb0 ≈ perturb0'
#
#	for p in eachslice(perturb, dims=3)
#
#		@show LinearAlgebra.norm(p)
#		
#		@show LinearAlgebra.norm(p-p') 
#
#	end 

#########################

	for dir1 in 1:2

		h = Bloch_Hamilt(w, dir1, P)

		data = WLO.get_wlo_data_mesh1(WLO.psiH_on_mesh1(nk, k0, perturb, h))
	
		set_results_onedir!(results, dir1, data)

	end 


	return results

end 



#===========================================================================#
#
# FoundFiles &  Read 
#
#---------------------------------------------------------------------------#


function FoundFiles(P::UODict; 
										target=nothing, get_fname::Function, kwargs...)::Bool

	FoundFiles0(get_fname(P), get_target(target; kwargs...))

end 



function FoundFiles0(Obs_fname::Function, 
										 observables::AbstractVector{<:AbstractString})::Bool

	!isempty(observables) && ReadWrite.FoundFiles_PhysObs(Obs_fname, observables, FILE_STORE_METHOD)

end





function Read(P::UODict; target=nothing, get_fname::Function, kwargs...)::Dict 

	Read0(get_fname(P), get_target(target; kwargs...))

end
	



function Read0(Obs_fname::Function, 
							 observables::AbstractVector{<:AbstractString}
												 )::Dict{String,Any}

	if isempty(observables)
		
		obs_found = [split(fn,".")[1] for fn in cd(readdir, Obs_fname())]

		return ReadWrite.Read_PhysObs(Obs_fname, obs_found, FILE_STORE_METHOD)

	else 

		return ReadWrite.Read_PhysObs(Obs_fname, observables, FILE_STORE_METHOD)

	end 

end 










































































































#############################################################################
end # module RibbonWLO 

