using Plots
using DelimitedFiles
using LaTeXStrings
using Printf
cd(@__DIR__)

include("../../plotting_params.jl")

figs_dir="../../../../figs/"

data_dir="../../../../data/supplimentary/crit_map_1024_1e8_bulcsu/"

method_params = Int.(readdlm(data_dir*"method_params.txt"))
T, Ttr, grid_size_tent,grid_size_crit = method_params
T_string::String = @sprintf "%.E" T
Ttr_string::String = @sprintf "%.E" Ttr
annotation_pos = (-0.25,  1.0)
markershapes = [:rect,:utriangle,:cross,:diamond,:circle]
line_alphas = reverse([1.0,0.7,0.5,0.3,0.1])
markersize = 3
msw = 0.2

#---------------------------critical map--------------------------

#-------------------dynamic renyi entropy spectrum--------------------
r=2.0
grid_size=grid_size_crit
#Os = Int.(readdlm(data_dir*"orders_crit.txt"))
Os = [1,2,4,8,10]
qs = 0.01:0.01:2
pl_crit_Kq = plot()
for (i,o) in enumerate(Os)
    files = readdir(data_dir)
    file_idx = findall(s -> occursin("critical_renyi-q_dq=0.01",s) && occursin("_order=$(o)",s),files)
    f_name = files[file_idx[1]]
    @show f_name, i, grid_size,r
    data = readdlm(data_dir * f_name)

    qs = data[:,1]
    Hs = data[:,2]
    plot!(pl_crit_Kq, qs, Hs, label=L"\tilde{K}_q(%$(o))", linewidth=1,markershape=:circle,mc=:red,ms=markersize,markerstrokecolor=:gray10,markerstrokewidth=msw,ma=line_alphas[i] ,la=line_alphas[i], color="red")

end

#=
#read crit longer time series, higher resolution
#data_crit_long = readdlm(data_dir2*"critical_renyi-q_dq=0.01_r=2.0_tmax=10^8_ttrans=10^6_grid=1024_order=10.dat")

#qs = data_crit_long[:,1]
#Hs = data_crit_long[:,2]
#plot!(pl_crit_Kq, qs[1:2:end], Hs[1:2:end], label=L"\tilde{K}_q(10), "*L"n=2^{10}, "*L"10^{8}", linewidth=3, alpha=1, color=:gray10)

#read crit, higher resolution
data_crit = readdlm(data_dir3*"critical_renyi-q_dq=0.01_r=2.0_tmax1E+07_ttrans1E+05_grid=1024_order=8.dat")

qs = data_crit[:,1]
Hs = data_crit[:,2]
plot!(pl_crit_Kq, qs, Hs, label=L"\tilde{K}_q(8), "*L"n=2^{10}, "*L"10^{7}",markershape=:circle,mc=:gray10,ms=3,markerstrokewidth=msw, linewidth=1, alpha=1, color=:gray10)
=#

plot!(pl_crit_Kq,[1.0],[0.5],st=:scatter,mc=:red,markerstrokewidth=1,markerstrokecolor=:gray10,ms=markersize,label="") #theoretical value at q=1


plot!(pl_crit_Kq, xlabel=L"q", ylabel=L"\tilde{K}_q(m)", ylim=[0,1.2], 
    xformatter=:auto,guidefontsize=guidefontsize, tickfontsize=tickfontsize-2, 
    legendfontsize=legendfontsize-2, dpi=300)
    #title = "Critical map ")
plot!(pl_crit_Kq, xlim=[0,2], ylim=[0.,1.1], yticks=[0,0.5,1.0], legend=:topright, left_margin=-1Plots.mm)
annotate!(pl_crit_Kq, annotation_pos, text("(d)", :left, 18))
#savefig(pl_crit_Kq, figs_dir*"critical_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.pdf")
#savefig(pl_crit_Kq, figs_dir*"critical_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.svg")


#-------------------static renyi entropy spectrum--------------------

#Os = Int.(readdlm(data_dir*"orders_tent_static.txt"))
Os = Int.(readdlm(data_dir*"orders_crit_static.txt"))[6:5:end-5] #don't show lower orders
qs = 0.01:0.01:2
pl_crit_Hq = plot()
for (i,o) in enumerate(Os)

    files = readdir(data_dir)
    file_idx = findall(s -> occursin("critical_static_renyi",s) && occursin("_order=$(o)",s),files)
    f_name = files[file_idx[1]]
    data = readdlm(data_dir * f_name)
    qs = data[:,1]
    Hs = data[:,2]
    plot!(pl_crit_Hq, qs, Hs/o, label=L"H_q(%$(o))/%$(o)", linewidth=1,markershape=:circle,mc=:red,ms=markersize,markerstrokewidth=msw, ma=line_alphas[i],la=line_alphas[i], color="red")
end
plot!(pl_crit_Hq, xlabel=L"q", ylabel=L"H_q(m)/m", 
    xformatter=:auto, guidefontsize=guidefontsize, tickfontsize=tickfontsize-2,
    legendfontsize=legendfontsize-2, dpi=300)
    #title = "Critical map ")

plot!(pl_crit_Hq,[1.0],[0.5],st=:scatter,mc=:red,markerstrokewidth=1,markerstrokecolor=:gray10,ms=markersize,label="") #theoretical value at q=1
plot!(pl_crit_Hq, xlim=[0,2], ylim=[0,1.1],yticks=[0,0.5,1.0], legend=:topright, left_margin=-1Plots.mm)
annotate!(pl_crit_Hq, annotation_pos, text("(c)", :left, 18))

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
grid_size=grid_size_tent
#Os = Int.(readdlm(data_dir*"orders_tent_dynamic.txt"))
Os = [1,2,4,8,10]
qs = 0.01:0.01:2
Hsa = analytical_renyi_entropy_spectrum(r, qs)
pl_tent_Kq = plot()
for (i,o) in enumerate(Os)

    files = readdir(data_dir)
    file_idx = findall(s -> occursin("asymmetric_triangular_renyi-q_dq=0.01",s) && occursin("_order=$(o)",s),files)
    f_name = files[file_idx[1]]
    data = readdlm(data_dir * f_name)


    qs = data[:,1]
    Hs = data[:,2]
    plot!(pl_tent_Kq, qs, Hs, label=L"\tilde{K}_q(%$(o))", linewidth=1,markershape=:circle,mc=:orange,ms=markersize,markerstrokewidth=msw, ma=line_alphas[i], la=line_alphas[i], color=:orange)
end
plot!(pl_tent_Kq, qs, Hsa, label=L"K_q", ls=:dash, lw=2, alpha=0.8, color="black")
plot!(pl_tent_Kq, xlabel="", ylabel=L"\tilde{K}_q(m), K_q", ylim=[0,1.2], 
    guidefontsize=guidefontsize, tickfontsize=tickfontsize-2, 
    legendfontsize=legendfontsize-2, dpi=300,xformatter=:none)
    #title = "Tent map " )
plot!(pl_tent_Kq, xlim=[0,2], ylim=[0.,1.1], yticks=[0,0.5,1], legend=:bottomleft, left_margin=-1Plots.mm)
annotate!(pl_tent_Kq, annotation_pos, text("(b)", :left, 18))
#savefig(pl_tent_Kq, figs_dir*"asymmetric_triangular_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.pdf")
#savefig(pl_tent_Kq, figs_dir*"asymmetric_triangular_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.svg")

#-------------------static renyi entropy spectrum--------------------

Os = Int.(readdlm(data_dir*"orders_tent_static.txt"))[6:5:end-5]
qs = 0.01:0.01:2
Hsa = analytical_renyi_entropy_spectrum(r, qs)
pl_tent_Hq = plot(dpi=300)
for (i,o) in enumerate(Os)

    files = readdir(data_dir)
    file_idx = findall(s -> occursin("asymmetric_triangular_static",s) && occursin("_order=$(o)",s),files)
    f_name = files[file_idx[1]]
    data = readdlm(data_dir * f_name)
    

    qs = data[:,1]
    Hs = data[:,2]
    plot!(pl_tent_Hq, qs, Hs/o, label=L"H_q(%$(o))/%$(o)", linewidth=1,markershape=:circle,mc=:orange,ms=markersize,markerstrokewidth=msw, ma=line_alphas[i], la=line_alphas[i], color=:orange)
end
plot!(pl_tent_Hq, qs, Hsa, label=L"K_q", ls=:dash, lw=2, alpha=0.8, color="black")
plot!(pl_tent_Hq, xlabel="", ylabel=L"H_q(m)/m", 
    guidefontsize=guidefontsize, tickfontsize=tickfontsize-2,
    legendfontsize=legendfontsize-2, dpi=300,xformatter=:none)
    #title = "Tent map ")
plot!(pl_tent_Hq, xlim=[0,2], ylim=[0,1.1], yticks=[0.0,0.5,1.0],legend=:bottomleft, left_margin=-1Plots.mm)
annotate!(pl_tent_Hq, annotation_pos, text("(a)", :left, 18))
#savefig(pl, figs_dir*"asymmetric_triangular_static_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.pdf")
#savefig(pl, figs_dir*"asymmetric_triangular_static_renyi_r=$(r)_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)_higherorder.svg")


pl = plot(pl_tent_Hq,pl_tent_Kq,pl_crit_Hq,pl_crit_Kq,layout=(2,2),size=(colfig_size[1],colfig_size[2]-300),left_margin=1Plots.mm,right_margin=2Plots.mm)
savefig(pl,figs_dir*"renyi_static_dynamic_entropy_spectra_tmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)"*"_1e8.png")