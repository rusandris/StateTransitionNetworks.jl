using DelimitedFiles
using Plots,LaTeXStrings
using Printf
include("plotting_params_main.jl")
include("plot_functions_chaotic_maps.jl")
#---------------------------------henon------------------------------------

#henon spectrum results
spectrums_dir = "data/henon_data/"

all_spec_files = readdir(spectrums_dir)
spectrum_files = all_spec_files[findall(f -> occursin("renyi_spectrums", f),all_spec_files)]
measure_files = all_spec_files[findall(f -> occursin("measures_q1", f),all_spec_files)]

all_spectrums = readdlm.(spectrums_dir .* spectrum_files) 
measures_q1 = readdlm.(spectrums_dir .* measure_files) 

param_file = all_spec_files[findall(f -> occursin("param_values", f),all_spec_files)][1]
special_ps = vec(readdlm(spectrums_dir * param_file))
orders = vec(readdlm(spectrums_dir * "orders.txt",Int64))


alphas = [1.0,1.0,1.0,1.0]

plot_params_spec = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
ylims=(-0.03,0.72),
titlefontsize=20,
left_margin=reduced_left_margin,
top_margin=top_margin,
legend= (0.7,0.9),
dpi=300,#size=(800,300)
xlabel=L"q",
yticks=[0.0,0.35,0.7],
yformatter=:none
)

qs::Vector{Float64} = all_spectrums[1][:,1]


#------------------henon ref values----------------------
ref_vals = [h_top_henon,λ_henon]


#plot only highest order spectrum
renyi_spectrums = [all_spectrums[i][:,5] for i in 1:length(all_spectrums)]
renyi_spectrums = hcat(renyi_spectrums...)
measures_q1_highest_order = [measures_q1[i][:,4][1] for i in 1:length(measures_q1)]

#labels = [L"\bar{K}_q(%$(orders[i]))" for i in 1:length(orders)]
labels = [L"a=%$(special_ps[i])" for i in 1:length(special_ps)]

pl_spec = plot_renyi_spectrums(qs,renyi_spectrums,measures_q1_highest_order,special_ps,ref_vals;
alphas=ones(4),
plot_ref_values=true,
labels = labels,
marker_shapes=marker_shapes,
marker_colors=marker_colors,
annotation_fontsize = annotation_fontsize,
annotation_factor = 1.1,
marker_size=marker_size_colored,
plot_params_spec...)

annotate!(pl_spec,subfigure_annotation_pos_henon, text("(b)", :left, 24))

#--------------------------measures fig---------------------
henon_results_dir = "data/henon_data/"
result_files = readdir(henon_results_dir)

ps_file = result_files[findfirst(f -> occursin("henon_p_values", f),result_files)]
ps = vec(readdlm(henon_results_dir * ps_file))

od_file = result_files[findfirst(f -> occursin("od", f),result_files)]
od_data = readdlm(henon_results_dir * od_file)
od_data = [od_data[i,:] for i in 1:size(od_data,1)]

orders_file = result_files[findfirst(f -> occursin("orders", f),result_files)]
orders = vec(readdlm(henon_results_dir * orders_file))


lyapexps_file = result_files[findfirst(f -> occursin("lyapexps", f),result_files)]
lyapexps_data = readdlm(henon_results_dir * lyapexps_file)

λs = lyapexps_data[:,2]

entropies_file = result_files[findall(f -> occursin("entropies", f),result_files)][1]
lambdas_file = result_files[findall(f -> occursin("lambdas", f),result_files)][1]

Ss = readdlm(henon_results_dir * entropies_file)[:,2:end]
Λs = readdlm(henon_results_dir * lambdas_file)[:,2:end]

#--------------plotting params------------
alphas = [0.2:0.2:1.0;]
#inset_indices = 1:201
inset_param_range = [1.2,1.25]
inset_ticks_entropies = [inset_param_range,[0.0,0.5]]
inset_ticks_lambdas = [inset_param_range,[0.0,0.5]]

s_labels = [L"S(%$(i))" for i in orders ]

Λ_labels = [L"\Lambda(%$(i))" for i in orders ]


plot_params_od = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
ylims=ylims_od,
yticks=[0.0,0.5,1.0],
xticks=[1.2,1.25,1.3,1.35,1.4],
titlefontsize=20,
left_margin=reduced_left_margin,
top_margin=reduced_top_margin,
#right_margin=4Plots.mm,
legend=:none,
#ylabel= L"x_n",
xformatter=:none,
yformatter=:none,
yguidefontrotation=0,
dpi=300)

plot_params_entropy = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
ylims=(-0.01,1.0),
yticks=[0.0,0.5,1.0],
xticks=[1.2,1.25,1.3,1.35,1.4],
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
ylims=(-0.01,1.2),
yticks=[0.0,0.5,1.0],
titlefontsize=20,
left_margin=reduced_left_margin,
#right_margin=4Plots.mm,
top_margin=reduced_top_margin,
legend=:none,
#ylabel= L"\Lambda",
xlabel = L"a",
vertical_lw=2,
xticks=[1.2,1.25,1.3,1.35,1.4],
yformatter=:none,
yguidefontrotation=0,
dpi=300)


pl_od_henon = plot_orbit_diagram(od_data,ps,special_ps;marker_shapes=marker_shapes,marker_offset=0.08,marker_size=marker_size_colored,marker_colors=marker_colors,plot_params_od...)
annotate!(pl_od_henon,subfigure_annotation_pos_henon, text("(d)", :left, 24))

pl_s_henon = plot_measure(ps,Ss,special_ps;
    labels = s_labels,
    λs = λs,
    alphas=alphas,
    inset_box=inset_box,
    inset_param_range=inset_param_range,
    marker_colors=marker_colors,
    orders=orders,
    inset_ticks=inset_ticks_entropies,
    inset_ylims = [-0.02,0.5],
    vertical_lw=2,
    marker_size=marker_size_colored,
    inset_tickfontsize=inset_tickfontsize,
    plot_params_entropy...)

annotate!(pl_s_henon[1],subfigure_annotation_pos_henon, text("(f)", :left, 24))

pl_lambda_henon = plot_measure(ps,Λs,special_ps;
    labels = Λ_labels,
    alphas=alphas,
    inset_box=inset_box,
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

annotate!(pl_lambda_henon[1],subfigure_annotation_pos_henon, text("(h)", :left, 24))

l = @layout [a{0.25h}; b{0.25h}; b{0.25h};d{0.25h}]
pl_henon = plot(pl_spec,pl_od_henon,pl_s_henon,pl_lambda_henon,layout = l,size=halfcolfig_size)
#savefig(pl_henon,"henon_4fig_spectrums_measures_T1e8_max_order8.png")
