module RibbonWLO  
#############################################################################

using OrderedCollections: OrderedDict 

import myLibs: ReadWrite, Utils, Lattices, TBmodel

import myLibs.Parameters: UODict 

import ..FILE_STORE_METHOD

import ..WLO, ..Helpers 

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
									 

calc_observables::Vector{String} = ["WannierBands1"]

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
#


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

function set_results_onedir!(results::AbstractDict, 
													dir1::Int,
													nus1::AbstractVector{<:Real},
													)::Nothing


	if haskey(results, "WannierBands1") 

		copy!(selectdim(results["WannierBands1"],2,dir1), nus1)

	end 

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

function Compute_(P::UODict, target, get_fname::Nothing=nothing; 
										kwargs...)::Dict{String,Any}

	nk = nr_kPoints(P) 

	w = width(P) 

	k0 = kPoint_start(P)

	results = init_results(w, get_target(target; kwargs...))

	@assert all_symms_preserved(P)

	for dir1 in 1:2 

		latt = lattice(w, 3-dir1)
	
		h = Matrix∘TBmodel.Bloch_Hamilt(Lattices.NearbyUCs(latt); nr_orb=4, Hopping=MODEL.get_hoppf(P))
	
		nus = WLO.get_wlo_data_mesh1(WLO.psiH_on_mesh1(nk, k0, h))

		set_results_onedir!(results, dir1, sort!(nus))

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

