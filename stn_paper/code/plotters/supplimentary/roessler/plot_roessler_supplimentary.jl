using DelimitedFiles
using Plots,LaTeXStrings
using Printf

cd(@__DIR__)
include("../../plotting_params.jl")
include("../../plot_functions_chaotic_maps.jl")

#----------------------------roessler----------------------------

#roessler results dir
roessler_results_dir = "../../../../data/supplimentary/roessler_data/"

#--------------------------measures fig---------------------

result_files = readdir(roessler_results_dir)

ps_file = result_files[findfirst(f -> occursin("roessler_p_values", f),result_files)]
ps = vec(readdlm(roessler_results_dir * ps_file))

od_file = result_files[findfirst(f -> occursin("od", f),result_files)]
od_data = readdlm(roessler_results_dir * od_file)
od_data = [od_data[i,:] for i in 1:size(od_data,1)]

orders_file = result_files[findfirst(f -> occursin("orders", f),result_files)]
orders = vec(readdlm(roessler_results_dir * orders_file))


lyapexps_file = result_files[findfirst(f -> occursin("lyapexps", f),result_files)]
lyapexps_data = readdlm(roessler_results_dir * lyapexps_file)

λs = lyapexps_data[:,2]

entropies_file = result_files[findall(f -> occursin("entropies", f),result_files)][1]
lambdas_file = result_files[findall(f -> occursin("lambdas", f),result_files)][1]

Ss = readdlm(roessler_results_dir * entropies_file)[:,2:end]
Λs = readdlm(roessler_results_dir * lambdas_file)[:,2:end]


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
xlims=[0.34,0.41],
yticks=[-9,-5.5,-2],
xticks=[0.35,0.36,0.37,0.38,0.39,0.4],
titlefontsize=20,
left_margin=reduced_left_margin,
#top_margin=reduced_top_margin,
#right_margin=4Plots.mm,
legend=:none,
#ylabel= L"x_n",
xformatter=:auto,
yformatter=:auto,
yguidefontrotation=0,
dpi=300)

plot_params_entropy = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
ylims=(-0.01,1.0),
yticks=[0.0,0.5,1.0],
xticks=[0.35,0.36,0.37,0.38,0.39,0.4],
titlefontsize=20,
left_margin=reduced_left_margin,
top_margin=reduced_top_margin,
#right_margin=4Plots.mm,
legend=:none,
#ylabel= L"S",
vertical_lw=2,
xformatter=:none,
yformatter=:none,
yguidefontrotation=0,
dpi=300)

plot_params_lambda = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
#framestyle=:box,
yticks=[0.0,0.5,1.0],
xticks=[0.35,0.36,0.37,0.38,0.39,0.4],
titlefontsize=20,
left_margin=reduced_left_margin,
right_margin=4Plots.mm,
top_margin=reduced_top_margin,
legend=:none,
#ylabel= L"\Lambda",
xlabel = L"a",
vertical_lw=2,
yformatter=:none,
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
    marker_colors=marker_colors,
    orders=orders,
    vertical_lw=2,
    marker_size=marker_size_colored,
    inset_tickfontsize=inset_tickfontsize,
    plot_params_entropy...)


pl_lambda_henon = plot_measure(ps,Λs,special_ps;
    labels = Λ_labels,
    alphas=alphas,
    inset_box=inset_box,
    inset_data = Λs_inset_data,
    inset_param_range=inset_param_range,
    marker_colors=marker_colors,
    orders=orders,
    inset_ticks=inset_ticks_entropies,
    inset_ylims = [-0.02,0.7],
    inset_tickfontsize=inset_tickfontsize,
    vertical_lw=2,
    marker_size=marker_size_colored,
    red_params=[1.22,1.227],
    plot_params_lambda...)

annotate!(pl_spec,subfigure_annotation_pos_two_col_inner, text("(b)", :left, annotation_fontsize))
annotate!(pl_od_henon,subfigure_annotation_pos_two_col_inner, text("(d)", :left, annotation_fontsize))
annotate!(pl_s_henon[1],subfigure_annotation_pos_two_col_inner, text("(f)", :left, annotation_fontsize))
annotate!(pl_lambda_henon[1],subfigure_annotation_pos_two_col_inner, text("(h)", :left, annotation_fontsize))


l = @layout [a{0.25h}; b{0.25h}; b{0.25h};d{0.25h}]
pl_henon = plot(pl_spec,pl_od_henon,pl_s_henon,pl_lambda_henon,layout = l)