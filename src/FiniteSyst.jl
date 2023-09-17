module FiniteSyst
##############################################################################


import myLibs.Parameters: UODict 

import myLibs: Lattices, BandStructure, Utils, ReadWrite, TBmodel, Operators

import myPlots 


import ..FILE_STORE_METHOD

import ..CalcWLO: MODEL 

import ..RibbonWLO



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


Dependencies = [MODEL]

usedkeys::Vector{Symbol} = [:width]


calc_observables::Vector{String} = ["kLabels", "Energy", "kLabel"]



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

linear_size = RibbonWLO.width 

function lattice(N::Int)::Lattices.Lattice 

	latt = Lattices.SquareLattice()

	Lattices.Superlattice!(latt, [N,N]; recenter=true)

	Lattices.ReduceDim!(latt) 
	
end 

lattice(P::UODict)::Lattices.Lattice = lattice(linear_size(P))

PosAtoms = Lattices.PosAtoms∘lattice 

#function Bloch_Hamilt(P::UODict)::Function 
#
#	latt = lattice(P)
#
#	Hopping = MODEL.get_hoppf(P)
#
#	l3 = Lattices.NearbyUCs(latt)
#
#	return TBmodel.Bloch_Hamilt(l3; nr_orb=4, Hopping=Hopping)
#
#end  

	

function get_BlochHamilt_argskwargs(P::UODict,
																		latt::Lattices.Lattice=lattice(P);
																		kwargs...
																		)
	
	hopp = MODEL.get_Hopping(P)

	for (k,v) in kwargs

		hopp[k] = v 

	end 

	return (Lattices.NearbyUCs(latt),), hopp 

end  


function get_BlochHamilt(args...; kwargs...)::Function 

	A, K = get_BlochHamilt_argskwargs(args...; kwargs...)

	return TBmodel.Bloch_Hamilt(A...; K...)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function get_operators(operator_names::AbstractVector{<:AbstractString},
											 atoms::AbstractMatrix{Float64};
											 nr_orb::Int, kw...
											)::Tuple{Vector{String}, Vector}

	kwargs = (nr_at=size(atoms,2), nr_orb=nr_orb, dim=2)

	operator_functions = map(operator_names) do n 

		d,f = myPlots.Transforms.parse_fstr_Cartesian(n)

		if d>0 

			return Operators.Position(d, atoms; kwargs..., fpos=f)

		end 

		for (N,S) in (("LocalPsi2",:LDOS),("IPR",:IPR))

			n==N && return getfield(Operators, S)(; kwargs...)

		end 

		error(n)

	end 

	return (operator_names, operator_functions)

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


#	obs = get_target(target; kwargs...) 

	obs = nothing 

	results = Compute_(P, obs; kwargs...)

	ReadWrite.Write_PhysObs(get_fname(P), FILE_STORE_METHOD, results)

	return results 
#	return isnothing(target) ? results : Utils.dict_keepkeys(results, obs) 

end 


function Compute_(P::UODict, target, get_fname::Nothing=nothing; 
									operators::AbstractVector{<:AbstractString},
										kwargs...)::Dict{String,Any}

	latt = lattice(P) 

	hopp = MODEL.get_Hopping(P)

	atoms = Lattices.PosAtoms(latt)


  return BandStructure.Diagonalize(

		TBmodel.Bloch_Hamilt(Lattices.NearbyUCs(latt,atoms); hopp...),
		
		zeros(1,1);
		
		dim=2,
		
		storemethod = FILE_STORE_METHOD,

		operators = get_operators(operators, atoms; hopp...),

		kwargs...
		
		)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#===========================================================================#
#
# FoundFiles &  Read 
#
#---------------------------------------------------------------------------#


function FoundFiles(P::UODict; 
										target=nothing, get_fname::Function, 
										operators::AbstractVector{<:AbstractString},
										kwargs...)::Bool

#	FoundFiles0(get_fname(P), 
#							calc_observables
##							get_target(target; kwargs...),
#							)

	ReadWrite.FoundFiles_PhysObs(get_fname(P), 
															 union(calc_observables, operators),
															 FILE_STORE_METHOD)

end





function Read(P::UODict; target=nothing, get_fname::Function, 
										operators::AbstractVector{<:AbstractString},
							kwargs...)::Dict 

#	Read0(get_fname(P), get_target(target; kwargs...))

	ReadWrite.Read_PhysObs(get_fname(P), 
												 union(calc_observables, operators),
												 FILE_STORE_METHOD)


end 

























































































































































































##############################################################################
end
