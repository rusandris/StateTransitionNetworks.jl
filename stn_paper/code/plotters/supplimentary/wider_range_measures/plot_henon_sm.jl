using DelimitedFiles
using Plots,LaTeXStrings
using Printf
using DynamicalSystems
#using ComplexityMeasures
cd(@__DIR__)
include("../../plotting_params.jl")
include("../../plot_functions_chaotic_maps.jl")

#--------------------------------------Henon--------------------------------------------

#read henon param values from spectrumdir
henon_results_dir = "../../../../data/supplimentary/henon_data_sm/"
results_files = readdir(henon_results_dir)
param_file = results_files[findfirst(f -> occursin("param_values", f),result_files)]
ps_henon = vec(readdlm(henon_results_dir * param_file))

global ps_wl_henon::Vector{Float64} = ps_henon[[1,3,4]] #[1.2265,1.27,1.4]
marker_shapes_wl = marker_shapes[[1,3,4]] #local plotting param


#---------------------------------walk length stats---------------------------

#----------------------------read in data-----------------------

result_files = readdir(henon_results_dir)

walk_lengths_file = result_files[findfirst(f -> occursin("walk_lengths", f),result_files)]
walk_lengths_henon = readdlm(henon_results_dir * walk_lengths_file)

walk_lengths_stats_file = result_files[findfirst(f -> occursin("wl_stats", f),result_files)]
means_henon,variances_henon,Ss_henon,Λs_henon = eachcol(readdlm(henon_results_dir * walk_lengths_stats_file))

pl_dist_henon = plot(xlabel=L"\mathcal{L}(t)",
#ylabel=L"p(\mathcal{L}(t))",
guidefontsize=guidefontsize,
tickfontsize=tickfontsize,
legendfontsize=legendfontsize,
legend=(0.2,0.9),
ylims=common_y_lims,
#yformatter=:none,
xticks = [0,200,400],
xlims = (-10,400),
dpi=300,
left_margin=reduced_left_margin,
right_margin=5Plots.mm,
yformatter=:none)

#L_ranges = [80:210,20:110,280:350]
L_ranges_henon = [50:270,170:230,300:380]

for (i,p) in enumerate(ps_wl_henon)
    histogram!(pl_dist_henon,walk_lengths_henon[:,i], 
    normalize=:pdf, 
    bins=nr_bins,
    color=wl_colors[i],
    alpha=histogram_alpha,
    lw = 0.00001,
    la=0.001,
    label="")

    μ = Ss_henon[i]*N_steps
    σ = sqrt(Λs_henon[i]*N_steps)

    #println("(S - <L>/t)/S =  ", (Ss_henon[i] - means_henon[i])/Ss_henon[i])
   # println("(Λ - Δ^2L/t)/S =  ", (Λs_henon[i] - variances_henon[i])/Λs_henon[i])

    plot!(pl_dist_henon, L_ranges_henon[i], gauss.(L_ranges_henon[i],μ,σ),
    color=wl_colors[i],
    linewidth=2,
    label="",
    dpi=300)

    #plot marker inside
    plot!([μ],[marker_pos_henon[i]],
    st=:scatter,
    ms=marker_size_colored,
    markerstrokewidth=0.5,
    mc=wl_colors[i],
    markershape=marker_shapes_wl[i],
    label=L"a=%$(p)")

end

#plot blue periodic
plot!(pl_dist_henon,[0.0,0.0],[0.0,0.6], lw=curve_lw,lc=:blue,label="")
plot!(pl_dist_henon,[0.0],[0.02],st=:scatter,
    mc=:blue,
    ms=marker_size_colored,
    markershape=:star5,
    markerstrokewidth=0.5,
    label=L"a=%$(ps_henon[2])")



#--------------------------measures fig---------------------
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

#--------------local plotting params------------

inset_box = bbox(0.2,0.0,0.3,0.25)
inset_param_range = [1.2,1.25]
inset_ticks_entropies = [inset_param_range,[0.0,0.5]]
inset_ticks_lambdas = [inset_param_range,[0.0,0.5]]

s_labels = [L"S(%$(i))" for i in orders ]

Λ_labels = [L"\Lambda(%$(i))" for i in orders ]


plot_params_od = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
ylims=(-0.5,1.7),
yticks=[-0.5,0.0,0.5,1.0,1.5],
xticks=[1.0,1.1,1.2,1.3,1.4],
titlefontsize=20,
left_margin=left_margin,
top_margin=reduced_top_margin_sm,
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
#framestyle=:box,
ylims=(-0.01,1.25),
yticks=[0.0,0.5,1.0],
xticks=[1.0,1.1,1.2,1.3,1.4],
titlefontsize=20,
left_margin=left_margin,
top_margin=reduced_top_margin_sm,
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
ylims=(-0.01,1.25),
yticks=[0.0,0.5,1.0],
titlefontsize=20,
left_margin=left_margin,
top_margin=reduced_top_margin_sm,
#right_margin=4Plots.mm,
legend=:none,
#ylabel= L"\Lambda",
xlabel = L"a",
vertical_lw=2,
xticks=[1.0,1.1,1.2,1.3,1.4],
yformatter=:none,
yguidefontrotation=0,
dpi=300)


pl_od_henon = plot_orbit_diagram(od_data,ps,ps_henon;ms_od = ms_od_henon,ma_od = ma_od_henon,
    marker_shapes=marker_shapes,marker_offset=0.17,marker_size=marker_size_colored,
    marker_colors=marker_colors,plot_params_od...)



pl_s_henon = plot_measure(ps,Ss,ps_henon;
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



pl_lambda_henon = plot_measure(ps,Λs,ps_henon;
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


annotate!(pl_dist_henon,subfigure_annotation_pos_bigfig_inner, text("(b)", :left, annotation_fontsize))
annotate!(pl_od_henon,subfigure_annotation_pos_bigfig_inner, text("(d)", :left, annotation_fontsize))
annotate!(pl_s_henon[1],subfigure_annotation_pos_bigfig_inner, text("(f)", :left, annotation_fontsize))
annotate!(pl_lambda_henon[1],subfigure_annotation_pos_bigfig_inner, text("(h)", :left, annotation_fontsize))

l = @layout [a{0.25h}; b{0.25h};c{0.25h};d{0.25h}]
pl_henon = plot(pl_dist_henon,pl_od_henon,pl_s_henon,pl_lambda_henon,layout = l)