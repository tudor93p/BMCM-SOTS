module MB  

import LinearAlgebra

import ..WLO 
import ..Helpers:Symmetries 
#import PyPlot 
import myLibs: Groups, Utils

import myLibs.Parameters: UODict

#import DelimitedFiles 

const quantized_wcc2_values = [0,0.5]

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function dispersion((kx,ky))

	1 - cos(kx) - cos(ky)

end 


function Delta((kx,ky),sigma)

	sigma*sin(kx) - im*sin(ky) 

end 

function singlet((s0,sx,sy),(kx,ky))

	s0 + sx*cos(kx) + sy*cos(ky)

end 

function complex_b((bx,by),sigma)

	bx + sigma*im*by 

end 

function H(k::AbstractVector{<:Real},
					 bs::AbstractVector{<:Real},
					 aux::Real=0
					 )::Matrix{ComplexF64}


	xi = dispersion(k) 

	Dp = Delta(k,1)
	
	Dm = Delta(k,-1)
	
	bp = complex_b(bs[1:2],1)
	
	bm = complex_b(bs[1:2],-1)

	sk = singlet(bs[3:5], k)


	return [[xi bm Dp sk];
				 [bp xi sk Dm];
				 [conj(Dp) sk -xi bp];
				 [sk conj(Dm) bm -xi]]

end 


function H(k::AbstractVector{<:Real},
					 bs::AbstractVector{<:Real},
					 perturb::AbstractMatrix{<:Number},
					 )::Matrix{ComplexF64}
	
	H(k, bs) .+= perturb 

end  


function H(k::AbstractVector{<:Real},
					 bs::AbstractVector{<:Real},
					 pstrength::Real,
					 perturb0::AbstractMatrix{<:Number},
					 )::Matrix{ComplexF64}

	LinearAlgebra.axpy!(pstrength, perturb0, H(k,bs))

end  


function get_s_param_from_b(b::Real,bp::Real)::Float64 
	
	max(0,abs(b)-abs(bp))*b

#	sqrt(max(0,abs(b)-abs(bp)))*b
end 



function bsss_cycle(theta::Real)::Vector{Float64}

	bx = -sin(theta)
	by = cos(theta)

	sx = get_s_param_from_b(by,bx)
	
	sy = -get_s_param_from_b(bx,by) 

	out = [bx,by,(sx+sy)*0.1,sx,sy]

	out .*= 0.3 
#	return 0.3*[bx,by,0,0,0]
#	return 0.3*[bx,by,0.01/0.3,sx,sy]
	return out

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

braiding_time(P::UODict)::Float64 = P[:braiding_time]  
braiding_theta(t::Real)::Float64  = t*2pi 
braiding_theta(P::UODict)::Float64 = braiding_theta(braiding_time(P))



function get_psiH(MBtime::Real, n::Int, k0::Real)::Array{ComplexF64,4}

	WLO.psiH_on_mesh(n, k0, H,  bsss_cycle(braiding_theta(MBtime)))

end 

function get_psiH(MBtime::Real, n::Int, k0::Real, 
									perturb0::AbstractArray{ComplexF64,4},
									strength::Real,
											)::Array{ComplexF64,4}

	iszero(strength) && return get_psiH(MBtime, n, k0)
	
	return WLO.psiH_on_mesh(n, k0, perturb0, 
													H,  bsss_cycle(braiding_theta(MBtime)),
													strength,
													)

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function get_Wx_Wy(th::Real,K::AbstractVector{<:Real}
									 )::NTuple{2,Matrix}
	
	Tuple(hcat((WLO.wlo1_evals_pm(k, d, H, bsss_cycle(th)) for k in K)...) for d=1:2)

end  
function get_Wxy_Wyx(th::Real,K::AbstractVector{<:Real}
										)::NTuple{2,Matrix}


	Tuple(hcat((WLO.wlo2_pm_evals([k,k], d, H, bsss_cycle(th))[:] for k in K)...) for d=[2,1])

end 


#function plot_Wx_Wy(ax, th::Real, K::AbstractVector{<:Real})
#
#	plot_Wx_Wy(ax, K, get_Wx_Wy(th, K))
#
#end 



function plot_Wx_Wy_old(ax, 
#									K::AbstractVector{<:Real}, Ws::NTuple{2,<:AbstractMatrix})
									th::Real, K::AbstractVector{<:Real})

######## dirty: alterantive data 
#
#	Ws = get_Wx_Wy(th, K)  
#
#
#	alt_path = "/media/tudor/Tudor/Work/2018_Higher-Order-Topology/MajoranaBraiding/Data/WLO1path/"
#
#	alt_Ws = map([0,1]) do d 
#		t = replace(string(round(th/2pi,digits=3)),'.'=>'p')
#		
#		t = rpad(t,5,'0') * "_$d"
#	        
#	#	)[:,:2].T
#	
#		t1 = joinpath(alt_path, t, "W_evals.dat")
#
#		while !isfile(t1)
#
#			println("Waiting for t='$t'")
#
#			sleep(10)
#
#		end 
#
#		@assert isfile(t1) t1 
#	
#		return DelimitedFiles.readdlm(t1)[:,1:2]
#	
#	end 
#
#
###########  


	for (xy,W) in enumerate(Ws)

		for pm=1:2, s=[-1,0,1] 

			alt_W = alt_Ws[xy][:,pm]

			alt_K = LinRange(extrema(K/pi)..., length(alt_W))


			kw = s==0 ? (label=string("\$","+-"[pm],"\$"),) : () 

			ax[xy].scatter(K/pi,W[pm,:].+s, s=0.4, c=["k","r"][pm],zorder=5,alpha=0.8; kw...)
			
			kw = (s==0 && pm==2) ? (label=string("P"),) : ()  
			
			ax[xy].scatter(alt_K, alt_W .+s, s=3.3, c="green", alpha=0.6, zorder=0; kw...)

		end 
		
		d1 = "xy"[xy]

		d2 = "xy"[3-xy] 


		ax[xy].set_xlabel("\$k_$d2/\\pi\$")

		ax[xy].set_ylabel("\$\\nu_$d1\$",rotation=0)

	@show size(W)

	error("dist_modulo,diff -> Utils.dist_periodic") 
	gap = max(1e-12, minimum(WLO.dist_modulo, diff(W,dims=1)))

		gap = gap<1e-10 ? "0.0" : string(round(gap, digits = Int(ceil(-log10(gap)))+2))[1:min(end,6)]

		gap = rpad(gap, 6, "0")

		ax[xy].set_title("\$W_$d1\$  (gap\$=$gap\$)")
		
		ax[xy].set_xlim(extrema(K/pi))

	end 

end 

function plot_Wxy_Wyx(ax, th::Real,K::AbstractVector{<:Real})


	for (xy,W) in enumerate(get_Wxy_Wyx(th, K))

		for pm=1:2, s=[-1,0,1]
		
			kw = s==0 ? (label=string("\$","+-"[pm],"\$"),) : () 

			ax[xy].scatter(K/pi, W[pm,:] .+s , c=["k","r"][pm], s=8; kw...)

		end 

		d1 = "xy"[xy]

		d2 = "xy"[3-xy] 

		ax[xy].set_xlabel("\$k_$d1/\\pi\$")

		ax[xy].set_ylabel("\$\\phi_{$d1$d2}\$",rotation=0) 
		
		ax[xy].set_title("\$W_{$d1$d2}\$",rotation=0) 

		ax[xy].set_xlim(extrema(K/pi))

	end 

end 



#function plot_Ws(theta::Real, filename::AbstractString)
#
#	PyPlot.ioff() 
#
#	fig,Ax = plot_Ws(theta)
#
#	fig.savefig(filename,format="png",dpi=300)
#
#	PyPlot.close()
#
#end 
#
#function plot_Wx_Wy(theta::Real, filename::AbstractString)
#
#	PyPlot.ioff() 
#
#	fig,Ax = plot_Wx_Wy(theta)
#
#	fig.savefig(filename,format="png",dpi=300)
#
#	PyPlot.close()
#
#end 
#
#function plot_Wx_Wy(theta::Real)
#
#	fig,Ax,K = init_fig_Ws(theta, [1,2]) 
#
#
#	plot_Wx_Wy(Ax, theta, K)
#	
#	for a in Ax 
#
#		a.legend(loc="upper right")
#
#		a.set_ylim(-.75,.75)
#
#	end 
#
#	fig.subplots_adjust(top=0.87,bottom=0.152,left=.09,right=0.98,
#											wspace=0.28)
#
#	return fig,Ax 
#
#end 
#
#
#
#function init_fig_Ws(theta::Real; nr_rc=[2,2], kwargs...)
#
#	figsize = (nr_rc[[2,1]]/maximum(nr_rc)) .* [1.2,1] * 6
#
#	fig,Ax = PyPlot.subplots(nr_rc...; figsize=figsize, kwargs...)
#
#
#	fig.suptitle(rpad(string("t = ", round(theta/2pi,digits=3)),9,'0') * " T")
#
#	fig.subplots_adjust(top=0.91,bottom=0.08,left=.09,right=0.98,
#										hspace=0.375,wspace=0.245)  
#
#	return fig, Ax, LinRange(-pi, pi, WLO.NR_kPOINTS*5)
#
#end 
#
#
#
#
#function plot_Ws(theta::Real)
#
#	fig, Ax, K = init_fig_Ws(theta) 
#
#	plot_Wx_Wy(Ax[1,:], theta, K)
#
#	plot_Wxy_Wyx(Ax[2,:], theta, K)
#
#
#	for a in Ax 
#
#		a.legend(loc="lower center")
#
#		a.set_ylim(-.7,.7)
#
#	end 
#
#
#	fig.subplots_adjust(top=0.91,bottom=0.08,left=.09,right=0.98,
#										hspace=0.375,wspace=0.245) 
#
#	return fig,Ax 
#
#end 
#
#

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function get_wcc1_data(args...; sectors=[+,-])
#
#	psiH = WLO.eigH_on_mesh(args...; halfspace=false)
#
#
#	W1xy, psiH = WLO.wlo1_on_mesh(1:2, Hparams...) 
#
#	W1eig = [[WLO.Wannier_subspace_on_mesh(W1,s) for s in sectors] for W1 in W1xy]
#
#	return W1eig,psiH
#
#end 
#
#
#function get_wcc1((W1eig,psiH), d1::Int)::Matrix{Float64}
#
#	vcat(reshape.(last.(W1eig[d1]),1,:)...)
#
#end 


#function get_wcc1gap(data)::Vector{Float64}
#
#	[get_wcc1gap(data, d1) for d1=1:2]
#
#end 
#
#
#function get_wcc1gap(data, d1::Int)::Float64
#
#	minimum(Utils.dist_periodic(eachrow(get_wcc1(data, d1))...,1))
#
#end 
#

#function get_wcc2((W1eig,psiH), d2::Int)::Matrix{Float64}
#
#	w2s = [WLO.wlo2_on_mesh(d2, w1wf, psiH) for (w1wf,) in W1eig[3-d2]]
#
#	w2 = cat(w2s..., dims=(1,2))
#
#
#	return reshape(WLO.store_on_mesh(WLO.get_periodic_eigvals, w2), size(w2,1),:)
#
#end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





MxRepr = kron(Groups.PauliMatrix(0),Groups.PauliMatrix(1))
MyRepr = kron(Groups.PauliMatrix(3),Groups.PauliMatrix(2))
CtRepr = kron(Groups.PauliMatrix(1),Groups.PauliMatrix(2))
PRepr = kron(Groups.PauliMatrix(1),Groups.PauliMatrix(3))
#MzRepr = kron  
TtRepr = kron(Groups.PauliMatrix(0),Groups.PauliMatrix(1))
TR2yRepr = kron(Groups.PauliMatrix(0),Groups.PauliMatrix(0))
TR2xRepr = kron(Groups.PauliMatrix(3),Groups.PauliMatrix(3))

fij_identity(k::NTuple{2,Int}, ::Function)::NTuple{2,Int} = k 



function operMx(h::AbstractMatrix)::Matrix

	Symmetries.UdAV(MxRepr,h)

end 

#function operMx((i,j)::T, n::Int, rev::Function)::Union{T,Nothing
#																							 } where T<:NTuple{2,Int}
#
#	i<=div(n,2)
#
#	(rev(i),j)
#
#end 



function operMx(k::AbstractVector)::Vector

	Symmetries.mirror(k, 1)

end 

function fij_operMx(k::NTuple{2,Int}, ind_minusk::Function
									 )::NTuple{2,Int} 
	
	Symmetries.mirror(k, 1, ind_minusk)

end 



function operMx((V,U)::Tuple{AbstractMatrix,AbstractMatrix},
								)::Matrix 

	Symmetries.UdAV(V,MxRepr,U)

end 


function operCt((V,U)::Tuple{AbstractMatrix,AbstractMatrix},
								)::Matrix 

	Symmetries(V,CtRepr,U)

end 

fij_operCt = fij_identity


function operMy(h::AbstractMatrix)::Matrix

	Symmetries.UdAV(MyRepr,h)

end 


function operMy(k::AbstractVector)::Vector

	Symmetries.mirror(k, 2)

end 


function fij_operMy(k::NTuple{2,Int}, ind_minusk::Function
									 )::NTuple{2,Int} 
	
	Symmetries.mirror(k, 2, ind_minusk)

end 




function operCt(h::AbstractMatrix)::Matrix

	-Symmetries.UdAV(CtRepr,h)

end 


function operCt(k::AbstractVector)::Vector

	copy(k)

end 



function operP(h::AbstractMatrix)::Matrix

	Symmetries(PRepr,-conj(h))

end 

fij_operP = Symmetries.inversion 


function operP(k::AbstractVector)::Vector

	Symmetries.inversion(k)

end  


function operTC2y(h::AbstractMatrix)::Matrix

	Symmetries.UdAV(TR2yRepr,conj(h))

end 


function operTC2x(h::AbstractMatrix)::Matrix

	Symmetries.UdAV(TR2xRepr, conj(h))

end 



function operTC2y(k::AbstractVector)::Vector

	out = Symmetries.inversion!(Symmetries.rotationC2(k,2))
	@assert out ==k .* [1,-1]
	return out 

end  


function fij_operTC2y(k::NTuple{2,Int}, ind_minusk::Function)::NTuple{2,Int}  

	Symmetries.inversion(Symmetries.rotationC2(k,2,ind_minusk), ind_minusk)

end 

function fij_operTC2x(k::NTuple{2,Int}, ind_minusk::Function)::NTuple{2,Int}  

	Symmetries.inversion(Symmetries.rotationC2(k,1,ind_minusk), ind_minusk)

end 


function operTC2x(k::AbstractVector)::Vector

	out = Symmetries.inversion!(Symmetries.rotationC2(k,1))
	@assert out ==k .* [-1,1]
	return out 

end 


function operTt(h::AbstractMatrix)::Matrix

	Symmetries.UdAV(TtRepr, conj(h))

end 


fij_operTt = Symmetries.inversion 


function operTt(k::AbstractVector)::Vector
	
	Symmetries.inversion(k)

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#






function symmetrize!(A::AbstractMatrix, 
										 op::Union{<:AbstractString, Symbol}, args...)::Nothing

	Symmetries.symmetrize!(A, getOpFun(op), args...)

end  

function sepSymmString(opers_::AbstractString, 
											 prefix::AbstractString="")::Vector{String}

	opers  = filter(!isempty, split(opers_,"+"))
	
	(isempty(opers) ||	in("None",opers)) && return String[]

	return [parse_operstring(op, prefix) for op in opers]

end  




function parse_operstring(op::AbstractString, 
											 prefix::AbstractString="")::String

	@assert !isempty(op) && !occursin("+",op)

	op[1:min(4,end)]=="oper" ? prefix * op : prefix * "oper" * op 
	
end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#catch e e isa UndefVarError || throw(e) 


function getOp_fij(op::AbstractString)::Function 
	
	fij1 = getOpFun(op, "fij_")

	function fij(n::Int, k0::Real)::Function 
		
		ind_minusk = WLO.ind_minusk(n, k0)

#	function fij(ind_minusk::Function)::Function 

		fij_(k::Vararg{Int,2})::NTuple{2,Int} = fij_(k) 

		fij_(k::NTuple{2,Int})::NTuple{2,Int} = fij1(k, ind_minusk)

	end 

	return fij 
	
end 

getOp_fij(op::AbstractString, args...)::Function = getOp_fij(op)(args...)


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


# dealing only with binary symmetries --- divide in halfspace 
# returns the indices {I} such that {I, f(I)} covers the entire space


function getOp_uniqueInds(op::AbstractString)::Function   

	f11 = getOpFun(op, "fij_")((1,1),-)

	return function uniqueInds_oper(n::Int, k0::Real
																	)::Tuple{<:AbstractVector{Int}, 
																					 <:AbstractVector{Int}
																					}
		ui0 = Base.OneTo(n-1)

		all(==(1),f11) && return (ui0, ui0)

		all(==(-1),f11) && return (WLO.uniqueInds_kMirror(n,k0), ui0)

		return Tuple(f==-1 ? WLO.uniqueInds_kMirror(n,k0) : ui0 for f=f11)

	end  

	
end 



function getOp_uniqueInds(op::AbstractString, n::Int, k0::Real
													 )::Tuple{<:AbstractVector{Int},
																		<:AbstractVector{Int}}
	
	getOp_uniqueInds(op)(n, k0)

end 



function getOpFun(op::AbstractString, args...)::Function 

	getOpFun(Symbol(parse_operstring(op, args...)))

end 


function getOpFun(op::Symbol)::Function 
	
	getproperty(@__MODULE__, op)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function symmetrize!(A::AbstractMatrix, opers::AbstractString,
										args...)::Nothing

	for op in sepSymmString(opers)
		
		symmetrize!(A, op, args...)

	end 

end 
	

function HamHasSymm(op::Union{AbstractString,Symbol},
										args...; kwargs...)::Bool 

	HamHasSymm(getOpFun(op), args...; kwargs...)

end 



function HamHasSymm(op::Function, fHamilt::Function, args...;
										nr_samples::Int64=10, kwargs...
									 )::Bool 

	all(eachcol(rand(2,nr_samples)*2pi)) do k 

		Symmetries.has_symm(fHamilt(op(k), args...), op, 
												fHamilt(k, args...);
												kwargs...)

	end 

end 

function has_symm(op::Union{AbstractString,Symbol},
									args...; kwargs...)::Bool

	Symmetries.has_symm(getOpFun(op), args...; kwargs...)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function calcplot_wcc12(theta::Real,
												kmesh::AbstractArray{<:Real};
												kwargs...)

	calcplot_wcc12(theta, kmesh, H; kwargs...)

end 


function calcplot_wcc12(theta::Real,
												kmesh::AbstractArray{<:Real},
												hamilt::Function,
												hargs...;
												fignum::Int=200,
												filename=nothing,
												savefig_kwargs=(),
												kwargs...
												)
	
	data = get_wcc1_data(hamilt, bsss_cycle(theta), hargs...)

	fig1, Ax1,  = init_fig_Ws(theta; num=fignum) 

	for (d1,d2) in [(1,2),(2,1)]

		WLO.plot_wcc(Ax1[1,d1], 
								 d1,
								 get_wcc1(data, d1),
								 view(selectdim(kmesh,1,d2),:);
								 kwargs...)

		WLO.plot_wcc2(Ax1[2,d1], d1,
									get_wcc2(data, d2),
									view(selectdim(kmesh,1,d1),:);
									kwargs...
									)

	end 

	isnothing(filename) || fig1.savefig(filename; savefig_kwargs...)


	return fig1,Ax1 

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




#
#
#function plot_y2stat(Y2Stat::AbstractArray{Float64,4},
#										 X::AbstractArray{<:Real},
#										 xaxis::Int;
#										 xlabel=nothing,
#										 suptitle=nothing,
#										 labels=nothing,
#										 kwargs...
#										 )
#
#	fig2,Ax2 = PyPlot.subplots(3,2; figsize=[1.2,1.5]*6, 
#														 sharex=true,sharey="row",
#														 kwargs...)
#
#	isnothing(suptitle) || fig2.suptitle(suptitle)
#
#	Ax2[1,1].set_ylabel("\$\\mu\$")
#	
#	Ax2[2,1].set_ylabel("\$\\sigma\$")
#
#	Ax2[3,1].set_ylabel("Wannier gap")
#
#
#	paxis = only(setdiff([1,4],xaxis))
#
#	paxis = Dict(4=>1,1=>3)[xaxis]
#	paxis2 = Dict(4=>1,1=>2)[xaxis]
#
#	
#	for (title,y2stat,ax) in zip(["Wxy","Wyx"],
#															 eachslice(Y2Stat,dims=3),
#															 eachcol(Ax2),
#															 )
#
#		for y in quantized_wcc2_values
#	
#			ax[1].plot(extrema(X),[y,y],c="gray",ls="--",lw=0.5, alpha=0.5)
#	
#		end 
#
#
#		for (y2s,lab) in zip(eachslice(y2stat,dims=paxis), 
#												 isnothing(labels) ? axes(y2stat,paxis) : labels)
#
#			for (a,Y) in zip(ax,eachslice(y2s,dims=paxis2))
#
#
#				a.plot(X, Y, label= isempty(lab) ? "None" : lab , alpha=0.7)
#
#			end 
#
#		end 
#	
#		ax[1].set_title(title)
#
#		isnothing(xlabel) || ax[end].set_xlabel(xlabel)
#
##		if size(y2stat,paxis)>1 
##
##			ax[1].legend()
##
##		end 
#
#	end  
#	
#		
#	Ax2[1,1].set_xlim(extrema(X))
#	Ax2[1,1].set_ylim(Utils.extend_limits(quantized_wcc2_values, 0.2))
#	
#	Ax2[2,1].set_ylim(0,0.04)
#
#	Ax2[3,1].set_ylim(0,0.08)
#	fig2.tight_layout() 
#
#
#	Ax2[1,1].legend()
#	Ax2[3,2].legend()
#
#	fig2.subplots_adjust(top=0.945, bottom=0.05, left=0.088, right=0.98, hspace=0.131, wspace=0.078)
#
#	return fig2,Ax2 
#
#
#end 






#############################################################################
end 
