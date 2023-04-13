import Dates  

t0 = Dates.now() 

import myLibs: ComputeTasks, Utils 

import BMCMSOTS  


include("input_file_10.jl")

tasks = [
				 init(BMCMSOTS,:CheckZero),
				 init(BMCMSOTS,:WannierBands1),
				 ];


ComputeTasks.missing_data.(tasks,show_missing=false);
error()

for t in tasks 

P, = ComputeTasks.get_first_paramcomb(t) 

for q in Base.product(((k=>v for v=input_checks[:allparams][k]) for k=[:kMesh_model, :preserved_symmetries])...)

	P1 = Utils.adapt_merge(P,q...,:nr_kPoints => 10)
	
	t.get_data(P1; force_comp=true, mute=false)

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


ComputeTasks.get_data_all.(tasks, 
													 shuffle=true, seed=4, 
													 mute=false,
#													 check_data=false,
													 )








