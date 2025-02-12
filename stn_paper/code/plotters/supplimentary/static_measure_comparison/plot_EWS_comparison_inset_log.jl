using DelimitedFiles
using Plots,LaTeXStrings
#cd("Documents/stn_research/STNResearch/Chaos/revision")
include("../../../StateTransitionNetworks.jl/stn_paper/code/plotters/plotting_params_main.jl")
include("../../../StateTransitionNetworks.jl/stn_paper/code/plotters/plot_functions_chaotic_maps.jl")
#---------------------------------load already computed results-------------------------

logistic_results_dir = "../prl_figure_data/logistic_1e8_reinit/"
result_files = readdir(logistic_results_dir)

ps_file = result_files[findfirst(f -> occursin("logistic_p_values", f),result_files)]
ps = vec(readdlm(logistic_results_dir * ps_file))

orders_file = result_files[findfirst(f -> occursin("orders", f),result_files)]
orders = vec(readdlm(logistic_results_dir * orders_file))

entropies_file = result_files[findall(f -> occursin("entrop", f),result_files)][1]
lambdas_file = result_files[findall(f -> occursin("lambda", f),result_files)][1]

spectrums_dir = "../prl_figure_data/spectrums/logistic1e8/"
special_ps = vec(readdlm(spectrums_dir * "param_values.txt"))


#---------------------------------od plot-----------------------------
#=
od_file = result_files[findfirst(f -> occursin("od", f),result_files)]
od_data = readdlm(logistic_results_dir * od_file)
od_data = [od_data[i,:] for i in 1:size(od_data,1)]
=#
recalc_dir = "SL_recalc_inset/"
recalc_files = readdir(recalc_dir)
od_file = recalc_files[findfirst(f -> occursin("logistic", f) && occursin("od", f) &&  occursin("nr_param_4001", f),recalc_files)]
od_data = readdlm(recalc_dir * od_file)
od_data = [od_data[i,:] for i in 1:size(od_data,1)]
ps_od = readdlm(recalc_dir*"params_od_logistic.txt")


linecolor=:gray10

plot_params_od = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
#framestyle=:box,
ylims=(0.0,1.45),
#xlims=(3.8,4.01),
yticks=[0.0,0.5,1.0],
titlefontsize=20,
left_margin=left_margin,
right_margin=reduced_right_margin ,
legend=:none,
ylabel= L"x_n",
yguidefontrotation=0,
xformatter=:none,
marker_size=marker_size_colored,
dpi=300)

pl_od_logistic = plot_orbit_diagram(od_data,ps_od,special_ps;ms_od = 0.3,ma_od=0.8,
    marker_shapes=marker_shapes,marker_offset=0.35,marker_size=marker_size_colored,
    marker_colors=marker_colors,plot_params_od...)

#---------------plot inset-------------------

inset_box_size = [0.2,0.4]
inset_box = bbox(0.65,0.05,inset_box_size...)
inset_xticks = [3.92,3.928]
inset_yticks = [0.4,1.0]

inset_xlims = [3.92,3.928]
inset_ylims = [0.4,1.0]
#plot!(pl_S, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,tickfontsize=inset_tickfontsize)
#plot!(pl_S[2],ps,od_data,ylims=inset_ylims,xlims=inset_xlims,lw=2,la=1.0,lc=linecolor,legend=false)
#lens!(pl_od_logistic,inset_xlims, inset_ylims, inset = (1,inset_box),grid=true,
#    framestyle=:box,lc=linecolor)
#plot!(pl_od_logistic[2],xticks=inset_xticks,yticks=inset_yticks)
#pl_od_logistic[2].series_list[:][:markersize].=1.8
#pl_od_logistic[2].series_list[:][:markercolor].=:red

#take highest order analytical results
Ss = readdlm(logistic_results_dir * entropies_file)[:,5]
Λs = readdlm(logistic_results_dir * lambdas_file)[:,5]

s_label = L"S(12)"
Λ_label = L"\Lambda(12)"

#--------------plotting params------------
alphas = [0.4:0.2:1.0;]
vertical_lw=2
linecolor =:gray10

plot_params = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
titlefontsize=20,
left_margin=left_margin+ 2Plots.mm,
right_margin=reduced_right_margin,
legend=:none,
xformatter=:none,
yguidefontrotation=0,
marker_size=marker_size_colored,
dpi=300)

#---------------------------S plot------------------------------
ylims = [0.0,1.0]
yticks = [0.0,0.5,1.0]
pl_S = plot(;ylabel= L"S",yticks=yticks,ylims=ylims,plot_params...)

#plot measure curves
plot!(pl_S,ps,Ss,lw=2,la=1.0,lc=linecolor,label = s_label)

#plot vertical lines on S plot
for (i,p) in enumerate(special_ps)
    plot!(pl_S,[p,p],[ylims[1],ylims[end]],ls=:dash,lw=vertical_lw,lc=marker_colors[i],label=false)
end

#---------------plot inset-------------------
ps_inset = readdlm("SL_recalc_inset/logistic_p_values_nr_param_161.txt")
S_inset = readdlm("SL_recalc_inset/logistic_entropies_T1E+08_Ttr1E+06_grid_32param_3.92_3.928.txt")[:,2]
inset_box = bbox(0.65,0.05,inset_box_size...)
inset_xticks = [3.92,3.928]
inset_yticks = [0.4,0.6]

inset_xlims = [3.92,3.928]
inset_ylims = [0.4,0.6]
#plot!(pl_S, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,tickfontsize=inset_tickfontsize)
#plot!(pl_S[2],ps_inset,S_inset,ylims=inset_ylims,xlims=inset_xlims,lw=2,la=1.0,lc=linecolor,legend=false)

lens!(pl_S,inset_xlims, inset_ylims, inset = (1,inset_box),grid=true,
    framestyle=:box,lc=linecolor)
plot!(pl_S[2],xticks=inset_xticks,yticks=inset_yticks)

#---------------------------Λ plot------------------------------
pl_Λ = plot(;ylabel= L"\Lambda",yticks=yticks,ylims=ylims,plot_params...)


#plot measure curves
plot!(pl_Λ,ps,Λs,lw=2,la=1.0,lc=linecolor,label = Λ_label)

#plot vertical lines on S plot
for (i,p) in enumerate(special_ps)
    plot!(pl_Λ,[p,p],[ylims[1],ylims[end]],ls=:dash,lw=vertical_lw,lc=marker_colors[i],label=false)
end

#---------------plot inset-------------------
Λ_inset = readdlm("SL_recalc_inset/logistic_lambda_T1E+08_Ttr1E+06_grid_32param_3.92_3.928.txt")[:,2]
inset_box = bbox(0.65,0.05,inset_box_size...)
inset_xticks = [3.92,3.928]
inset_yticks = [0.02,0.22]

inset_xlims = [3.92,3.928]
inset_ylims = [0.02,0.22]
#plot!(pl_Λ, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,tickfontsize=inset_tickfontsize)
#plot!(pl_Λ[2],ps_inset,Λ_inset,ylims=inset_ylims,xlims=inset_xlims,lw=2,la=1.0,lc=linecolor,legend=false)

lens!(pl_Λ,inset_xlims, inset_ylims, inset = (1,inset_box),grid=true,
    framestyle=:box,lc=linecolor)
plot!(pl_Λ[2],xticks=inset_xticks,yticks=inset_yticks)

EWS_standard_log = readdlm("var_acf_log_T1e8_TTr_1e6_rs_3.8_4.txt")
#-------------------------var plot---------------------------
ylims=[0,0.2]
yticks=[0.0,0.1,0.2]
var = EWS_standard_log[:,1]
pl_var = plot(;ylabel=L"\sigma^2",yticks=yticks,ylims=ylims,plot_params...)

#plot measure curves
plot!(pl_var,ps,var,lw=2,la=1.0,lc=linecolor,label = L"\sigma^2")

#plot vertical lines on var plot
for (i,p) in enumerate(special_ps)
    plot!(pl_var,[p,p],[ylims[1],ylims[end]],ls=:dash,lw=vertical_lw,lc=marker_colors[i],label=false)
end

#---------------plot inset-------------------
EWS_standard_log_inset = readdlm("var_acf_log_zoom_T1e8_TTr_1e6_rs_3.92_3.928.txt")
var_inset = EWS_standard_log_inset[:,2]
inset_box = bbox(0.65,0.5,inset_box_size...)
inset_xticks = [3.92,3.928]
inset_yticks = [0.08,0.12]

inset_xlims = [3.92,3.928]
inset_ylims = [0.08,0.12]
#plot!(pl_var, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,tickfontsize=inset_tickfontsize)
#plot!(pl_var[2],ps_inset,var_inset,ylims=inset_ylims,xlims=inset_xlims,lw=2,la=1.0,lc=linecolor,legend=false)

lens!(pl_var,inset_xlims, inset_ylims, inset = (1,inset_box),grid=true,
    framestyle=:box,lc=linecolor)
plot!(pl_var[2],xticks=inset_xticks,yticks=inset_yticks)


#-------------------------ac plot---------------------------
ylims=[-1.0,0.1]
yticks=[-1.0,-0.5,0.1]
ac = EWS_standard_log[:,2]
pl_ac = plot(;ylabel=L"ACF(1)",ylims=ylims,yticks=yticks,plot_params...,xformatter=:auto)


#plot measure curves
plot!(pl_ac,ps,ac,lw=2,la=1.0,lc=linecolor,label = L"ACF(1)")

#plot vertical lines on var plot
for (i,p) in enumerate(special_ps)
    plot!(pl_ac,[p,p],[ylims[1],ylims[end]],ls=:dash,lw=vertical_lw,lc=marker_colors[i],label=false)
end

#---------------plot inset-------------------
ac_inset = EWS_standard_log_inset[:,3]
inset_box = bbox(0.65,0.5,inset_box_size...)
inset_xticks = [3.92,3.928]
inset_yticks = [-0.5,-0.3]

inset_xlims = [3.92,3.928]
inset_ylims = [-0.5,-0.3]
#plot!(pl_ac, inset=inset_box, subplot=2, framestyle=:box,xticks=inset_xticks,yticks=inset_yticks,tickfontsize=inset_tickfontsize)
#plot!(pl_ac[2],ps_inset,ac_inset,ylims=inset_ylims,xlims=inset_xlims,lw=2,la=1.0,lc=linecolor,legend=false)

lens!(pl_ac,inset_xlims, inset_ylims, inset = (1,inset_box),grid=true,
    framestyle=:box,lc=linecolor)
plot!(pl_ac[2],xticks=inset_xticks,yticks=inset_yticks)

#l = @layout [a{0.25h}; b{0.25h}; b{0.25h};d{0.25h}]
l = (5,1)
pl_logistic = plot(pl_od_logistic,pl_S,pl_Λ,pl_var,pl_ac,layout = l,size=(600,1000))
#savefig(pl_henon,"EWS_comparison_logistic.png")