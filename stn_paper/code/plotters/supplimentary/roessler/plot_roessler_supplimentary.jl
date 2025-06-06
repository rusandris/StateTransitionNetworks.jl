using DelimitedFiles
using Plots,LaTeXStrings
using Printf

cd(@__DIR__)
include("../../plotting_params.jl")
include("../../plot_functions_chaotic_maps.jl")

#----------------------------roessler----------------------------

#roessler results dir
roessler_results_dir = "../../../../data/supplimentary/roessler_data/"

fig_dir_name = "figs"
fig_dir = "../../../../" * fig_dir_name * "/" 
mkpath(fig_dir)
method_params = readdlm(roessler_results_dir*"method_params.txt")
T,Ttr,grid_size = Int.(method_params[1:3])
T_string::String = @sprintf "%.E" T
Ttr_string::String = @sprintf "%.E" Ttr
grid_edges = Int.(method_params[4,:])
#--------------------------measures fig---------------------

result_files = readdir(roessler_results_dir)

ps_file = result_files[findfirst(f -> occursin("roessler_p_values", f),result_files)]
ps = vec(readdlm(roessler_results_dir * ps_file))

od_file = result_files[findfirst(f -> occursin("od_roessler", f),result_files)]
od_data = readdlm(roessler_results_dir * od_file)
od_data = [od_data[i,:] for i in 1:size(od_data,1)]

orders_file = result_files[findfirst(f -> occursin("orders", f),result_files)]
orders = vec(readdlm(roessler_results_dir * orders_file))

avg_cross_time_file = result_files[findfirst(f -> occursin("avg", f),result_files)]
avg_cross_times = readdlm(roessler_results_dir * avg_cross_time_file)

lyapexps_file = result_files[findfirst(f -> occursin("lyapexps", f),result_files)]
lyapexps_data = readdlm(roessler_results_dir * lyapexps_file)

λs = lyapexps_data[:,2] .* avg_cross_times[:,2]  #.* 5.49

entropies_file = result_files[findall(f -> occursin("entropies", f),result_files)][1]
lambdas_file = result_files[findall(f -> occursin("lambdas", f),result_files)][1]

Ss = readdlm(roessler_results_dir * entropies_file)[:,2:end] #[end:-1:1,2:end] #revert columns for plotting
Λs = readdlm(roessler_results_dir * lambdas_file)[:,2:end] #revert columns for plotting


#--------------plotting params------------
alphas = [0.2:0.2:1.0;]
#inset_indices = 1:201

s_labels = [L"S(%$(i))" for i in orders ]
Λ_labels = [L"\Lambda(%$(i))" for i in orders ]


plot_params_od = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
ylims=(-9,-2),
#xlims=[0.34,0.41],
yticks=[-8,-5.0,-2],
xticks=[0.35,0.36,0.37,0.38,0.39,0.4],
titlefontsize=20,
left_margin=reduced_left_margin,
#top_margin=reduced_top_margin,
#right_margin=4Plots.mm,
legend=:none,
xaxis=:flip,
title="Rössler system, PSOS " * L"y=2",
ylabel= L"x_{\tau}",
xformatter=:none,
yformatter=:auto,
yguidefontrotation=0,
dpi=300)

plot_params_entropy = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
ylims=(-0.1,1.0),
yticks=[0.0,0.5,1.0],
xticks=[0.35,0.36,0.37,0.38,0.39,0.4],
titlefontsize=20,
left_margin=reduced_left_margin,
top_margin=reduced_top_margin,
#right_margin=4Plots.mm,
legend=:topright,
xaxis=:flip,
ylabel= L"S,\langle T ~ \rangle \lambda",
vertical_lw=2,
xformatter=:none,
#yformatter=:none,
yguidefontrotation=0,
dpi=300)

plot_params_lambda = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
#framestyle=:box,
ylims = (-0.1,2.1),
yticks=[0.0,1.0,2.0,3.0],
xticks=[0.35,0.36,0.37,0.38,0.39,0.4],
titlefontsize=20,
left_margin=reduced_left_margin,
right_margin=4Plots.mm,
#top_margin=reduced_top_margin,
#legend=:none,
legend=:topright,
ylabel= L"\Lambda",
xlabel = L"b",
xaxis=:flip,
vertical_lw=2,
#yformatter=:none,
yguidefontrotation=0,
dpi=300)

special_ps = [0.36,0.368,0.4] #dummy 

pl_od_roessler = plot_orbit_diagram(od_data,ps,special_ps;ms_od = ms_od_henon,ma_od = ma_od_henon,
    marker_shapes=marker_shapes_roessler,marker_offset=0.3,marker_size=marker_size_colored,
    marker_colors=marker_colors_roessler,plot_params_od...)


pl_s_roessler = plot_measure(ps,Ss,special_ps;
    labels = s_labels,
    λs = λs,
    alphas=alphas,
    inset=false,
    marker_colors=marker_colors_roessler,
    orders=orders,
    vertical_lw=2,
    marker_size=marker_size_colored,
    inset_tickfontsize=inset_tickfontsize,
    plot_params_entropy...)


pl_lambda_roessler = plot_measure(ps,Λs,special_ps;
    labels = Λ_labels,
    alphas=alphas,
    inset=false,
    marker_colors=marker_colors_roessler,
    orders=orders,
    vertical_lw=2,
    marker_size=marker_size_colored,
    plot_params_lambda...)

annotate!(pl_od_roessler,(1.13,1.0), text("(a)", :left, annotation_fontsize))
annotate!(pl_s_roessler,(1.13,1.0), text("(b)", :left, annotation_fontsize))
annotate!(pl_lambda_roessler,(1.13,1.0), text("(c)", :left, annotation_fontsize))

l = @layout [a{0.33h}; b{0.33h}; c{0.33h}]
pl_roessler = plot(pl_od_roessler,pl_s_roessler,pl_lambda_roessler,layout = l,size=colfig_size)

savefig(pl_roessler,fig_dir * "roessler_supplimentary_Ttmax$T_string"*"_ttrans$Ttr_string"*"_grid=$(grid_size)"*".png")
