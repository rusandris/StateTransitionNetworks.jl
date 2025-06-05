using Plots
using DelimitedFiles
using LaTeXStrings
using Printf
cd(@__DIR__)

include("../../plotting_params.jl")

figs_dir="../../../../figs/"

data_dir="../../../../data/supplimentary/crit_tent_maps/"

method_params = Int.(readdlm(data_dir*"method_params.txt"))
T, Ttr, grid_size = method_params
T_string::String = @sprintf "%.E" T
Ttr_string::String = @sprintf "%.E" Ttr
#---------------------------critical map--------------------------

#-------------------dynamic renyi entropy spectrum--------------------
r=2.0

Os = Int.(readdlm(data_dir*"orders_crit.txt"))
qs = 0.01:0.01:2
pl_crit_Kq = plot()
for (i,o) in enumerate(Os)
    f_name = data_dir*"critical_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(o).dat"
    @show f_name, i, grid_size,r
    data = readdlm(f_name)
    qs = data[:,1]
    Hs = data[:,2]
    plot!(pl_crit_Kq, qs, Hs, label=L"\tilde{K}_q(%$(o))", linewidth=3, alpha=(1/(length(Os)+1))*i, color="red")
end
plot!(pl_crit_Kq, xlabel="", ylabel=L"\tilde{K}_q(m)", ylim=[0,1.2], 
    xformatter=:none,guidefontsize=guidefontsize, tickfontsize=tickfontsize, 
    legendfontsize=legendfontsize-2, dpi=300,
    title = "Critical map " * L"(r=2,10^7, n = 2^5)")
plot!(pl_crit_Kq, xlim=[0,2], ylim=[0.,4.], yticks=[0:1:4;], legend=:topright, left_margin=3Plots.mm)
#annotate!(pl, (-0.18,  1), text("(a)", :left, 24))
#savefig(pl_crit_Kq, figs_dir*"critical_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.pdf")
#savefig(pl_crit_Kq, figs_dir*"critical_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.svg")


#-------------------static renyi entropy spectrum--------------------

#Os = Int.(readdlm(data_dir*"orders_tent_static.txt"))
Os = Int.(readdlm(data_dir*"orders_crit.txt")) 
qs = 0.01:0.01:2
pl_crit_Hq = plot()
for (i,o) in enumerate(Os)
    f_name = data_dir*"critical_static_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(o).dat"
    @show f_name, i, grid_size
    data = readdlm(f_name)
    qs = data[:,1]
    Hs = data[:,2]
    plot!(pl_crit_Hq, qs, Hs/o, label=L"H_q(%$(o))/%$(o)", linewidth=3, alpha=(1/(length(Os)+1))*i, color="red")
end
plot!(pl_crit_Hq, xlabel="", ylabel=L"H_q(m)/m", 
    xformatter=:none, guidefontsize=guidefontsize, tickfontsize=tickfontsize,
    legendfontsize=legendfontsize-2, dpi=300,
    title = "Critical map " * L"(r=2,10^7, n = 2^5)")
plot!(pl_crit_Hq, xlim=[0,2], ylim=[0,4], legend=:topright, left_margin=3Plots.mm)
#plot!(pl, ylim=[0,2.5])
#savefig(pl, figs_dir*"critical_static_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid$(grid_size)_higherorder.pdf")
#savefig(pl, figs_dir*"critical_static_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.svg")


#----------------------------tent map----------------------------------

function analytical_renyi_entropy_spectrum(w, qs; verbose=true)
    Hs = zeros(length(qs))
    for (i,q) in enumerate(qs)
        verbose && @show q
        Hs[i] = analytical_renyi_entropy(w, q)
    end
    return Hs
end

function analytical_renyi_entropy(w, q)
    if q==0
        return log(2)
    elseif q==1
        return -w*log(w)-(1-w)*log(1-w)
    else
        return log(w^q+(1-w)^q)/(1-q)
    end
end

#-------------------dynamic renyi entropy spectrum--------------------
r=0.8

Os = Int.(readdlm(data_dir*"orders_tent_dynamic.txt"))
qs = 0.01:0.01:2
Hsa = analytical_renyi_entropy_spectrum(r, qs)
pl_tent_Kq = plot()
for (i,o) in enumerate(Os)
    f_name = data_dir*"asymmetric_triangular_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(o).dat"
    @show f_name, i, grid_size
    data = readdlm(f_name)
    qs = data[:,1]
    Hs = data[:,2]
    plot!(pl_tent_Kq, qs, Hs, label=L"\tilde{K}_q(%$(o))", linewidth=3, alpha=(1/(length(Os)+1))*i, color=:orange)
end
plot!(pl_tent_Kq, qs, Hsa, label=L"K_q", ls=:dash, lw=2, alpha=0.8, color="black")
plot!(pl_tent_Kq, xlabel=L"q", ylabel=L"\tilde{K}_q(m), K_q", ylim=[0,1.2], 
    guidefontsize=guidefontsize, tickfontsize=tickfontsize, 
    legendfontsize=legendfontsize-2, dpi=300,
    title = "Tent map " * L"(r=0.8,10^7, n = 2^5)")
plot!(pl_tent_Kq, xlim=[0,2], ylim=[0.,1.1], yticks=[0,0.5,1], legend=:bottomleft, left_margin=3Plots.mm)
#annotate!(pl, (-0.18,  1), text("(a)", :left, 24))
#savefig(pl_tent_Kq, figs_dir*"asymmetric_triangular_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.pdf")
#savefig(pl_tent_Kq, figs_dir*"asymmetric_triangular_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.svg")

#-------------------static renyi entropy spectrum--------------------


Os = Int.(readdlm(data_dir*"orders_tent_static.txt"))
qs = 0.01:0.01:2
Hsa = analytical_renyi_entropy_spectrum(r, qs)
pl_tent_Hq = plot(dpi=300)
for (i,o) in enumerate(Os)
    f_name = data_dir*"asymmetric_triangular_static_renyi-q_dq=0.01_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_order=$(o).dat"
    @show f_name, i, grid_size
    data = readdlm(f_name)
    qs = data[:,1]
    Hs = data[:,2]
    plot!(pl_tent_Hq, qs, Hs/o, label=L"H_q(%$(o))/%$(o)", linewidth=3, alpha=(1/(length(Os)+1))*i, color=:orange)
end
plot!(pl_tent_Hq, qs, Hsa, label=L"K_q", ls=:dash, lw=2, alpha=0.8, color="black")
plot!(pl_tent_Hq, xlabel=L"q", ylabel=L"H_q(m)/m", 
    guidefontsize=guidefontsize, tickfontsize=tickfontsize,
    legendfontsize=legendfontsize-2, dpi=300,
    title = "Tent map " * L"(r=0.8,10^7, n = 2^5)")
plot!(pl_tent_Hq, xlim=[0,2], ylim=[0,1.1], yticks=[0.0,0.5,1.0],legend=:bottomleft, left_margin=3Plots.mm)
#annotate!(pl, (-0.18,  1), text("(a)", :left, 24))
#savefig(pl, figs_dir*"asymmetric_triangular_static_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.pdf")
#savefig(pl, figs_dir*"asymmetric_triangular_static_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.svg")


pl = plot(pl_crit_Hq,pl_crit_Kq,pl_tent_Hq,pl_tent_Kq,layout=(2,2),size=(colfig_size[1],colfig_size[2]-300),left_margin=10Plots.mm,right_margin=10Plots.mm)
savefig(pl,figs_dir*"renyi_static_dynamic_entropy_spectra_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)"*".png")