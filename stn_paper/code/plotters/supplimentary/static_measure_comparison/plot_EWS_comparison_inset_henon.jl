using DelimitedFiles
using Plots,LaTeXStrings
#cd("Documents/stn_research/STNResearch/Chaos/revision")
include("../../../StateTransitionNetworks.jl/stn_paper/code/plotters/plotting_params_main.jl")
include("../../../StateTransitionNetworks.jl/stn_paper/code/plotters/plot_functions_chaotic_maps.jl")
#---------------------------------load already computed results-------------------------

henon_results_dir = "../prl_figure_data/henon1e8/"
result_files = readdir(henon_results_dir)

ps_file = result_files[findfirst(f -> occursin("henon_p_values", f),result_files)]
ps = vec(readdlm(henon_results_dir * ps_file))

orders_file = result_files[findfirst(f -> occursin("orders", f),result_files)]
orders = vec(readdlm(henon_results_dir * orders_file))

entropies_file = result_files[findall(f -> occursin("entrop", f),result_files)][1]
lambdas_file = result_files[findall(f -> occursin("lambda", f),result_files)][1]

spectrums_dir = "../prl_figure_data/spectrums/henon1e8/"
special_ps = vec(readdlm(spectrums_dir * "henon_param_values.txt"))



#---------------------------------od plot-----------------------------
#=
od_file = result_files[findfirst(f -> occursin("od", f),result_files)]
od_data = readdlm(henon_results_dir * od_file)
od_data = [od_data[i,:] for i in 1:size(od_data,1)]
=#

recalc_dir = "SL_recalc_inset/"
recalc_files = readdir(recalc_dir)
od_file = recalc_files[findfirst(f -> occursin("henon", f) && occursin("od", f) &&  occursin("nr_param_4001", f),recalc_files)]
od_data = readdlm(recalc_dir * od_file)
od_data = [od_data[i,:] for i in 1:size(od_data,1)]
ps_od = readdlm(recalc_dir*"params_od_henon.txt")

plot_params_od = (
guidefontsize=guidefontsize+5,
legendfontsize=legendfontsize+5,
tickfontsize=tickfontsize+2,
#framestyle=:box,
#ylims=(-1.4,1.6),
#yticks=[-1.0,0.0,1.0],
yticks=[0.7,1.0,1.3],
xticks=[1.2,1.3,1.4],
titlefontsize=20,
left_margin=left_margin,
bottommargin=10Plots.mm,
#right_margin=4Plots.mm,
legend=:none,
ylabel= L"x_n",
#xlabel= L"a",
xformatter=:none,
#yformatter=:none,
ylims=(0.7,1.4),
size=(1000,300),
yguidefontrotation=0,
dpi=300)

pl_od_henon = plot_orbit_diagram(od_data,ps_od,special_ps;ms_od = 0.5,ma_od=0.9,
    marker_shapes=marker_shapes,marker_offset=0.06,marker_size=marker_size_colored,
    marker_colors=marker_colors,plot_params_od...)

#---------------plot inset-------------------
inset_box_size = [0.2,0.4]
inset_box = bbox(0.65,0.05,inset_box_size...)
inset_xticks = [1.305,1.32]
inset_yticks = [1.2,1.4]

inset_xlims = [1.305,1.32]
inset_ylims = [1.2,1.4]
#plot!(pl_S, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,tickfontsize=inset_tickfontsize)
#plot!(pl_S[2],ps,Ss,ylims=inset_ylims,xlims=inset_xlims,lw=2,la=1.0,lc=linecolor,legend=false)

#lens!(pl_od_henon,inset_xlims, inset_ylims, inset = (1,inset_box),grid=true,
#    framestyle=:box,lc=linecolor)#tickfontsize=inset_tickfontsize)
#plot!(pl_od_henon[2],xticks=inset_xticks,yticks=inset_yticks)


#take highest order analytical results
Ss = readdlm(henon_results_dir * entropies_file)[:,5]
Λs = readdlm(henon_results_dir * lambdas_file)[:,5]

s_label = L"S(12)"
Λ_label = L"\Lambda(12)"

#--------------plotting params------------
alphas = [0.4:0.2:1.0;]
inset_box = bbox(0.65,0.13,inset_box_size...)
vertical_lw=2
linecolor =:gray10

plot_params = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
titlefontsize=20,
left_margin=reduced_left_margin-4Plots.mm,
right_margin=reduced_right_margin,
legend=:none,
xformatter=:none,
yguidefontrotation=0,
marker_size=marker_size_colored,
dpi=300)

#---------------------------S plot------------------------------
ylims = [0.0,1.0]
yticks = [0.0,0.5,1.0]
pl_S = plot(;yticks=yticks,ylims=ylims,plot_params...)
#ylabel= L"S"
#plot measure curves
plot!(pl_S,ps,Ss,lw=2,la=1.0,lc=linecolor,label = s_label)

#plot vertical lines on S plot
for (i,p) in enumerate(special_ps)
    plot!(pl_S,[p,p],[ylims[1],ylims[end]],ls=:dash,lw=vertical_lw,lc=marker_colors[i],label=false)
end

#---------------plot inset-------------------
ps_inset = readdlm("SL_recalc_inset/henon_p_values_nr_param_301.txt")
S_inset = readdlm("SL_recalc_inset/henon_entropies_T1E+08_Ttr1E+06_b_0.3_grid_32param_1.305_1.32.txt")[:,2]
inset_box = bbox(0.65,0.1,inset_box_size...)
inset_xticks = [1.305,1.32]
inset_yticks = [0.2,0.4]

inset_xlims = [1.305,1.32]
inset_ylims = [0.2,0.4]
#plot!(pl_S, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,tickfontsize=inset_tickfontsize)
#plot!(pl_S[2],ps_inset,S_inset,ylims=inset_ylims,xlims=inset_xlims,lw=2,la=1.0,lc=linecolor,legend=false)

lens!(pl_S,inset_xlims, inset_ylims, inset = (1,inset_box),grid=true,
    framestyle=:box,lc=linecolor)#tickfontsize=inset_tickfontsize)
plot!(pl_S[2],xticks=inset_xticks,yticks=inset_yticks)


#---------------------------Λ plot------------------------------
pl_Λ = plot(;yticks=yticks,ylims=ylims,plot_params...)
#ylabel= L"\Lambda"

#plot measure curves
plot!(pl_Λ,ps,Λs,lw=2,la=1.0,lc=linecolor,label = Λ_label)

#plot vertical lines on S plot
for (i,p) in enumerate(special_ps)
    plot!(pl_Λ,[p,p],[ylims[1],ylims[end]],ls=:dash,lw=vertical_lw,lc=marker_colors[i],label=false)
end

#---------------plot inset-------------------
Λ_inset = readdlm("SL_recalc_inset/henon_lambdas_T1E+08_Ttr1E+06_b_0.3_grid_32param_1.305_1.32.txt")[:,2]
inset_box = bbox(0.65,0.1,inset_box_size...)
inset_xticks = [1.305,1.32]
inset_yticks = [0.1,1.2]

inset_xlims = [1.305,1.32]
inset_ylims = [0.1,1.2]
#plot!(pl_Λ, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,tickfontsize=inset_tickfontsize)
#plot!(pl_Λ[2],ps_inset,Λ_inset,ylims=inset_ylims,xlims=inset_xlims,lw=2,la=1.0,lc=linecolor,legend=false)
lens!(pl_Λ,inset_xlims, inset_ylims, inset = (1,inset_box),grid=true,
    framestyle=:box,lc=linecolor)
plot!(pl_Λ[2],xticks=inset_xticks,yticks=inset_yticks)

EWS_standard_henon = readdlm("var_acf_henon_T1e8_TTr_1e6_rs_3.8_4.txt")
#-------------------------var plot---------------------------
ylims=[0.43,0.6]
yticks=[0.43,0.5,0.6]
var = EWS_standard_henon[:,1]
pl_var = plot(;yticks=yticks,ylims=ylims,plot_params...)
#ylabel=L"\sigma^2"

#plot measure curves
plot!(pl_var,ps,var,lw=2,la=1.0,lc=linecolor,label = L"\sigma^2")

#plot vertical lines on var plot
for (i,p) in enumerate(special_ps)
    plot!(pl_var,[p,p],[ylims[1],ylims[end]],ls=:dash,lw=vertical_lw,lc=marker_colors[i],label=false)
end

#---------------plot inset-------------------
EWS_standard_henon_inset = readdlm("var_acf_henon_zoom_T1e8_TTr_1e6_as_1.305_1.32_4.txt")
var_inset = EWS_standard_henon_inset[:,2]
inset_box = bbox(0.65,0.1,inset_box_size...)
inset_xticks = [1.305,1.32]
inset_yticks = [0.45,0.55]

inset_xlims = [1.305,1.32]
inset_ylims = [0.45,0.55]
#plot!(pl_var, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,tickfontsize=inset_tickfontsize)
#plot!(pl_var[2],ps_inset,var_inset,ylims=inset_ylims,xlims=inset_xlims,lw=2,la=1.0,lc=linecolor,legend=false)
lens!(pl_var,inset_xlims, inset_ylims, inset = (1,inset_box),grid=true,
    framestyle=:box,lc=linecolor)
plot!(pl_var[2],xticks=inset_xticks,yticks=inset_yticks)
#-------------------------ac plot---------------------------
ylims=[-0.7,0.5]
yticks=[-0.7,-0.1,0.5]
ac = EWS_standard_henon[:,3]
pl_ac = plot(;ylims=ylims,yticks=yticks,plot_params...,xformatter=:auto)
#ylabel=L"ACF(1)"

#plot measure curves
plot!(pl_ac,ps,ac,lw=2,la=1.0,lc=linecolor,label = L"ACF(1)")

#plot vertical lines on var plot
for (i,p) in enumerate(special_ps)
    plot!(pl_ac,[p,p],[ylims[1],ylims[end]],ls=:dash,lw=vertical_lw,lc=marker_colors[i],label=false)
end

#---------------plot inset-------------------
ac_inset = EWS_standard_henon_inset[:,4]
inset_box = bbox(0.65,0.1,inset_box_size...)
inset_xticks = [1.305,1.32]
inset_yticks = [-0.5,-0.35]

inset_xlims = [1.305,1.32]
inset_ylims = [-0.5,-0.35]
#plot!(pl_ac, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,tickfontsize=inset_tickfontsize)
#plot!(pl_ac[2],ps_inset,ac_inset,ylims=inset_ylims,xlims=inset_xlims,lw=2,la=1.0,lc=linecolor,legend=false)

lens!(pl_ac,inset_xlims, inset_ylims, inset = (1,inset_box),grid=true,
    framestyle=:box,lc=linecolor)
plot!(pl_ac[2],xticks=inset_xticks,yticks=inset_yticks)

#l = @layout [a{0.25h}; b{0.25h}; b{0.25h};d{0.25h}]
l = (5,1)
pl_henon = plot(pl_od_henon,pl_S,pl_Λ,pl_var,pl_ac,layout = l,size=(600,1000))
#savefig(pl_henon,"EWS_comparison_henon.png")