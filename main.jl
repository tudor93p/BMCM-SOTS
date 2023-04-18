import Dates  

t0 = Dates.now() 

using Distributed 

import myLibs: ComputeTasks, Utils 

import BMCMSOTS  


include("input_file_32.jl")

tasks = [
				 init(BMCMSOTS,:CheckZero),
#				 init(BMCMSOTS,:WannierBands1),
				 ];
	
ComputeTasks.missing_data.(tasks,show_missing=false)#gethostname()=="tudor-HP")



#@assert !in(gethostname(),["tudor-HP","horon"])


prep_all = ComputeTasks.get_data_all_prep.(tasks, 
#													 shuffle=true, 
													 seed=4, 
													 mute=false,
#													 check_data=false,
													 )

function relevant_params(criterion::Function, vals::AbstractVector) 

	length(vals)==1 && return vals 

	A = map(criterion,vals)

	return vals[filter!(!isnothing,[findfirst(A),findfirst(!,A)])]

end   

function relevant_params(::Val{:preserved_symmetries})

	relevant_params(==("All"),
									input_checks[:allparams][:preserved_symmetries] 
									)
end  

relevant_params(::Val{:nr_kPoints})=filter(>(3),sort(unique([5, nworkers()])))


function relevant_params(::Val{:kMesh_model})
	
	relevant_params(==("uniform")∘lowercase,
									get(input_checks[:allparams],:kMesh_model,["Uniform"])
									)
end 

relevant_params(k::Symbol)=[k=>v for v=relevant_params(Val(k))]

#relevant_params(k::Vector{Symbol}) = map(relevant_params, k)

relevant_params(k::Symbol...) = map(relevant_params, k)




for (t,(active,a,k)) in zip(tasks,prep_all)

	active || continue 	

	P, = ComputeTasks.get_first_paramcomb(t) 

	for q in Base.product(relevant_params(:kMesh_model, :preserved_symmetries, :nr_kPoints)...)

		t.get_data(Utils.adapt_merge(P,q...); force_comp=true, mute=false)
	
	end 



end 


#ComputeTasks.get_data_one.(tasks, mute=false);

#shuffle = gethostname()=="tudor-HP"



#@info "Preparations finished. Proceed to calculations? y/n" 

#if occursin("y",lowercase(readline(stdin)))  


t1 = Dates.DateTime("2023-02-27T18:11:46.870")

if Dates.now()<t1 

	println("Calculations start at: ",t1)

	while Dates.now() < t1 
	
		sleep(0.1)
	
	end 

end 


for (t,(active,args,kwargs)) in zip(tasks,prep_all)

	ComputeTasks.get_data_all(t, active, args...; kwargs...)

end 








