using DelimitedFiles
using Plots,LaTeXStrings
using Printf
cd(@__DIR__)
include("../plotting_params.jl")
include("../plot_functions_chaotic_maps.jl")
#---------------------------------logistic------------------------------------

#-----------------------------logistic spectrum results------------------------
logistic_results_dir = "../../../data/main/logistic_data/"

result_files = readdir(logistic_results_dir)
spectrum_files = result_files[findall(f -> occursin("renyi_spectrums", f),result_files)]
measure_files = result_files[findall(f -> occursin("measures_q1", f),result_files)]

all_spectrums = readdlm.(logistic_results_dir .* spectrum_files) 
measures_q1 = readdlm.(logistic_results_dir .* measure_files) 

param_file = result_files[findall(f -> occursin("param_values", f),result_files)][1]
special_ps = vec(readdlm(logistic_results_dir * param_file))
orders = vec(readdlm(logistic_results_dir * "orders.txt",Int64))

#marker_colors = fill(:red,4)
#alphas = [0.2,0.4,0.6,1.0]
#alphas = reverse(alphas)
alphas = [1.0,1.0,1.0,1.0]

plot_params_spec = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
ylims=(-0.03,0.72),
yticks=[0.0,0.35,0.7],
titlefontsize=20,
left_margin=left_margin,
top_margin=top_margin,
legend= (0.7,0.8),
dpi=300,
xlabel=L"q",
ylabel=L"\tilde{K}_q,K_q",
yguidefontrotation=0,
)

qs::Vector{Float64} = all_spectrums[1][:,1]

#------------------logistic ref values----------------------

ref_vals = [log(2),log(2)]


#plot only highest order spectrum
renyi_spectrums = [all_spectrums[i][:,5] for i in 1:length(all_spectrums)]
renyi_spectrums = hcat(renyi_spectrums...)
measures_q1_highest_order = [measures_q1[i][:,4][1] for i in 1:length(measures_q1)]

#labels = [L"\bar{K}_q(%$(orders[i]))" for i in 1:length(orders)]
labels = [L"r=%$(special_ps[i])" for i in 1:length(special_ps)]

pl_spec = plot_renyi_spectrums(qs,renyi_spectrums,measures_q1_highest_order,special_ps,ref_vals;
alphas=ones(4),
plot_ref_values=false,
labels = labels,
marker_shapes=marker_shapes,
marker_colors=marker_colors,
annotation_fontsize = annotation_fontsize,
annotation_factor = 0.9,
plot_params_spec...)



#plot analytical renyi for logistic map r=4 (log2)
plot!(pl_spec,[qs[1],qs[end]],[log(2),log(2)],ls=:dot,lw=2,lc=:gray10,label=L"K_q")

#--------------------------measures fig---------------------

result_files = readdir(logistic_results_dir)

ps_file = result_files[findfirst(f -> occursin("logistic_p_values", f),result_files)]
ps = vec(readdlm(logistic_results_dir * ps_file))

od_file = result_files[findfirst(f -> occursin("od", f),result_files)]
od_data = readdlm(logistic_results_dir * od_file)
od_data = [od_data[i,:] for i in 1:size(od_data,1)]

orders_file = result_files[findfirst(f -> occursin("orders", f),result_files)]
orders = vec(readdlm(logistic_results_dir * orders_file))

lyapexps_file = result_files[findfirst(f -> occursin("lyap", f),result_files)]
lyapexps_data = readdlm(logistic_results_dir * lyapexps_file)

λs = lyapexps_data

entropies_file = result_files[findall(f -> occursin("entrop", f),result_files)][1]
lambdas_file = result_files[findall(f -> occursin("lambda", f),result_files)][1]

Ss = readdlm(logistic_results_dir * entropies_file)[:,2:end]
Λs = readdlm(logistic_results_dir * lambdas_file)[:,2:end]



#--------------plotting params------------
alphas = [0.4:0.2:1.0;]
#inset_indices = 301:351
inset_param_range = [3.8,3.88]
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
xticks=[3.8,3.85,3.9,3.95,4.0],
titlefontsize=20,
top_margin=reduced_top_margin,
left_margin=left_margin,
right_margin=reduced_right_margin ,
legend=:none,
ylabel= L"x_n",
yguidefontrotation=0,
xformatter=:none,
marker_size=marker_size_colored,
dpi=300)

plot_params_entropy = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
#framestyle=:box,
ylims=(-0.01,1.0),
#xlims=(3.8,4.01),
yticks=[0.0,0.5,1.0],
xticks=[3.8,3.85,3.9,3.95,4.0],
titlefontsize=20,
left_margin=left_margin,
right_margin=reduced_right_margin,
top_margin=reduced_top_margin,
legend=:none,
ylabel= L"S,\lambda",
xformatter=:none,
yguidefontrotation=0,
marker_size=marker_size_colored,
dpi=300)

plot_params_lambda = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
#framestyle=:box,
ylims=(-0.01,1.25),
#xlims=(3.8,4.01),
yticks=[0.0,0.5,1.0],
titlefontsize=20,
left_margin=left_margin,
right_margin=reduced_right_margin,
top_margin=reduced_top_margin,
legend=:none,
ylabel= L"\Lambda",
xlabel = L"r",
xticks=[3.8,3.85,3.9,3.95,4.0],
yguidefontrotation=0,
marker_size=marker_size_colored,
dpi=300)


pl_od_logistic = plot_orbit_diagram(od_data,ps,special_ps;ms_od = ms_od_logistic,ma_od = ma_od_logistic,
    marker_shapes=marker_shapes,marker_offset=0.35,
    marker_size=marker_size_colored,marker_colors=marker_colors,plot_params_od...)



pl_s_logistic = plot_measure(ps,Ss,special_ps;
    labels = s_labels,
    λs = λs,
    alphas=alphas,
    inset_box=inset_box,
    inset_param_range=inset_param_range,
    marker_shapes=marker_shapes,
    marker_colors=marker_colors,
    orders=orders,
    inset_ticks=inset_ticks_entropies,
    inset_ylims = [-0.02,0.5],
    inset_tickfontsize = inset_tickfontsize,
    vertical_lw=2,
    marker_size=marker_size_colored,
    plot_params_entropy...)



pl_lambda_logistic = plot_measure(ps,Λs,special_ps;
    labels = Λ_labels,
    alphas=alphas,
    inset_box=inset_box,
    #inset_rw_vals=Λs_rw,
    inset_param_range=inset_param_range,
    marker_shapes=marker_shapes,
    marker_colors=marker_colors,
    orders=orders,
    inset_ticks=inset_ticks_entropies,
    inset_ylims = [-0.02,0.7],
    vertical_lw=2,
    marker_size=marker_size_colored,
    inset_tickfontsize=inset_tickfontsize,
    red_params=[3.82,3.83],
    plot_params_lambda...)

annotate!(pl_spec,subfigure_annotation_pos_two_col, text("(a)", :left, annotation_fontsize))
annotate!(pl_od_logistic,subfigure_annotation_pos_two_col, text("(c)", :left, annotation_fontsize))
annotate!(pl_s_logistic[1],subfigure_annotation_pos_two_col, text("(e)", :left, annotation_fontsize))
annotate!(pl_lambda_logistic[1],subfigure_annotation_pos_two_col, text("(g)", :left, annotation_fontsize))

l = @layout [a{0.25h}; b{0.25h}; b{0.25h};d{0.25h}]
pl_logistic = plot(pl_spec,pl_od_logistic,pl_s_logistic,pl_lambda_logistic,layout = l)
#savefig(pl_logistic,"logistic_4fig_spectrums_measures_test_T1e8_max_order12.png")
