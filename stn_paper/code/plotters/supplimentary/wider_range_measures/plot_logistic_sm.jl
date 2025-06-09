using DelimitedFiles
using Plots,LaTeXStrings
using Printf
using DynamicalSystems
#using ComplexityMeasures
cd(@__DIR__)
include("../../plotting_params.jl")
include("../../plot_functions_chaotic_maps.jl")

#---------------------------logistic------------------------------------

#-----------------------walk length stats---------------------------

#------------------------read in data-----------------------------


logistic_results_dir = "../../../../data/supplimentary/logistic_data_sm/"
result_files = readdir(logistic_results_dir)

param_file = result_files[findfirst(f -> occursin("param_values", f),result_files)]

global ps_log::Vector{Float64} = vec(readdlm(logistic_results_dir * param_file))
global ps_wl_log::Vector{Float64} = ps_log[[1,3,4]] #[3.82842,3.85,4.0]

walk_lengths_file = result_files[findfirst(f -> occursin("walk_lengths", f),result_files)]
walk_lengths_log = readdlm(logistic_results_dir * walk_lengths_file)

walk_lengths_stats_file = result_files[findfirst(f -> occursin("wl_stats", f),result_files)]
means_log,variances_log,Ss_log,Λs_log = eachcol(readdlm(logistic_results_dir * walk_lengths_stats_file))

marker_shapes_wl = marker_shapes[[1,3,4]] #local plotting param

pl_dist_log = plot(xlabel=L"\mathcal{L}(t)",
ylabel=L"p(\mathcal{L}(t))",
guidefontsize=guidefontsize,
tickfontsize=tickfontsize,
legendfontsize=legendfontsize,
legend=(0.5,0.9),
ylims=common_y_lims,
yticks = [0.0,0.03,0.06],
xticks = [0,250,500],
xlims = (-10,500),
dpi=300,
left_margin=left_margin)

L_ranges_log = [0:200,140:180,440:500]

for (i,p) in enumerate(ps_wl_log)

    histogram!(pl_dist_log,walk_lengths_log[:,i], 
    normalize=:pdf, 
    bins=nr_bins,
    color=wl_colors[i],
    alpha=histogram_alpha,
    label="",
    lw = 0.00001,
    la=0.001)

    μ = Ss_log[i]*N_steps
    σ = sqrt(Λs_log[i]*N_steps)

    #println("(S - <L>/t)/S =  ", (Ss_log[i] - means_log[i])/Ss_log[i])
    #println("(Λ - Δ^2L/t)/S =  ", (Λs_log[i] - variances_log[i])/Λs_log[i])

    plot!(pl_dist_log, L_ranges_log[i], gauss.(L_ranges_log[i],μ,σ),
    color=wl_colors[i],
    linewidth=curve_lw,
    label="",
    dpi=300)

    #plot marker inside
    plot!(pl_dist_log,[μ],[marker_pos_log[i]],
    st=:scatter,
    ms=marker_size_colored,
    markerstrokewidth=0.5,
    mc=wl_colors[i],
    markershape=marker_shapes_wl[i],
    label=L"r=%$(p)")

end

#plot blue periodic
plot!(pl_dist_log,[0.0,0.0],[0.0,0.6],lc=:blue, lw=curve_lw,label="")
plot!(pl_dist_log,[0.0],[0.02],st=:scatter,
    mc=:blue,
    ms=marker_size_colored,
    markershape=:star5,
    markerstrokewidth=0.5,
    label=L"r=%$(ps_log[2])")


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


entropies_file = result_files[findall(f -> occursin("entro", f),result_files)][1]
lambdas_file = result_files[findall(f -> occursin("lambda", f),result_files)][1]

Ss = readdlm(logistic_results_dir * entropies_file)[:,2:end]
Λs = readdlm(logistic_results_dir * lambdas_file)[:,2:end]

#--------------local plotting params------------
inset_box = bbox(0.3,0.0,0.3,0.25)
inset_param_range = [3.8,3.88]
inset_ticks_entropies = [inset_param_range,[0.0,0.5]]
inset_ticks_lambdas = [inset_param_range,[0.0,0.5]]

s_labels = [L"S(%$(i))" for i in orders ]

Λ_labels = [L"\Lambda(%$(i))" for i in orders ]

plot_params_od = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
#framestyle=:box,
ylims=(-0.5,1.7),
#xlims=(3.8,4.01),
yticks=[-0.5,0.0,0.5,1.0,1.5],
titlefontsize=20,
left_margin=left_margin,
right_margin=reduced_right_margin,
top_margin=reduced_top_margin_sm,
legend=:none,
ylabel= L"x_\tau",
yguidefontrotation=0,
xformatter=:none,
xticks=[3.5,3.625,3.75,3.875,4.0],
marker_size=marker_size_colored,
dpi=300)

plot_params_entropy = (
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
top_margin=reduced_top_margin_sm,
legend=:none,
ylabel= L"S,\lambda",
xticks=[3.5,3.625,3.75,3.875,4.0],
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
top_margin=reduced_top_margin_sm,
legend=:none,
ylabel= L"\Lambda",
xlabel = L"r",
xticks=[3.5,3.625,3.75,3.875,4.0],
yguidefontrotation=0,
marker_size=marker_size_colored,
dpi=300)


pl_od_logistic = plot_orbit_diagram(od_data,ps,ps_log;ms_od = ms_od_logistic,ma_od = ma_od_logistic,
    marker_shapes=marker_shapes,marker_size=marker_size_colored,
    marker_colors=marker_colors,plot_params_od...)


pl_s_logistic = plot_measure(ps,Ss,ps_log;
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



pl_lambda_logistic = plot_measure(ps,Λs,ps_log;
    labels = Λ_labels,
    alphas=alphas,
    inset_box=inset_box,
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


annotate!(pl_dist_log,subfigure_annotation_pos_one_col, text("(a)", :left, annotation_fontsize))
annotate!(pl_od_logistic,subfigure_annotation_pos_one_col, text("(c)", :left, annotation_fontsize))
annotate!(pl_s_logistic[1],subfigure_annotation_pos_one_col, text("(e)", :left, annotation_fontsize))
annotate!(pl_lambda_logistic[1],subfigure_annotation_pos_one_col, text("(g)", :left, annotation_fontsize))


l = @layout [a{0.25h}; b{0.25h};c{0.25h};d{0.25h}]
pl_logistic = plot(pl_dist_log,pl_od_logistic,pl_s_logistic,pl_lambda_logistic,layout = l)
#savefig(pl_logistic,"logistic_4fig_wl_measures_test_T1e8_max_order12.png")