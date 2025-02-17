using DelimitedFiles
using CSV,DataFrames
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting_params.jl")


#find annotations that mark change of signal type
#fibrillation onsets can fall between two RR annotations
#put fib onset between two discrete timepoints (integer indices)?
#put fib onset on the i-th RR interval which begins before the fib annotation
function find_annotations(annotation,sig_change_data,times)

    ann_idx = 1 #indx of annotation in annotation list

    ann_times = Float64[]
    rr_idxs = Int64[]

    #find every specified annotation (onset or termination) flag
    while true 
        ann_idx = findnext(s-> s==annotation,sig_change_data[:,end],ann_idx+1)
        isnothing(ann_idx) && break

        ann_time = sig_change_data[ann_idx,2]/128 #convert to time 
        push!(ann_times,ann_time)

        rr_idx = findlast(t -> t <= ann_time,times)
        push!(rr_idxs,rr_idx)

    end
    return rr_idxs,ann_times
end  

#--------------------read dir----------------------
data_dir = "../../../../data/supplimentary/fibrillation_ECG/fibrillation_results/"
result_files = readdir(data_dir)
samples = readdlm(data_dir*"sample_ids.txt",String)
method_params = readdlm(data_dir*"method_params.txt")
grid_size,w,window_size,grid_min,grid_max = method_params
grid_size,w,window_size = Int.([grid_size,w,window_size])
grid_edges = (grid_min,grid_max)
rrs_dir = "../../../../data/supplimentary/fibrillation_ECG/Long_Term_AF_Database/"
#--------------------plot data----------------------

plot_params = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
titlefontsize=20,
left_margin=reduced_left_margin,
top_margin=reduced_top_margin,
right_margin=reduced_right_margin,
legend=:none,
xformatter=:none,
yformatter=:none,
yguidefontrotation=0,
marker_size=marker_size_colored,
dpi=300)

#plots are zoomed in on onset
pre_onset_offset = 500 #xaxis limits: nr of idxs before onset
post_onset_offset = 200 #xaxis limits: nr of idxs before onset
ann_onset = "(AFIB"
ann_term = "(N"

xticks = [[1600:300:2500;],[4900:300:6000;],[1200:300:2300;],[1600:300:2800;]]

plots_grid = []
plots_OP = []
for (i,sample) in enumerate(samples)
    pl_rr = plot(;legend=false,ylims=grid_edges,yticks=[grid_edges[1]:0.4:grid_edges[2];],plot_params...)
    pl_S = plot(;legend=false,ylims=(-0.1,2.5),yticks=[0.0,1.0,2.0],plot_params...)
    pl_L = plot(;legend=false,ylims=(-0.1,3.0),yticks=[0:1:3;],plot_params...,) 
    pl_var = plot(;legend=false,ylims=(-0.02,0.08),yticks=[0.0:0.03:0.06;],plot_params...)
    pl_ac = plot(;legend=false,ylims=(-0.8,1.2),yticks=[-0.5,0.0,0.5,1.0],plot_params...,xformatter=:auto)
    xlabel!(pl_ac,"index")

    if i == 1
        plot!(pl_rr,ylabel=L"RR",yformatter=:auto)
        plot!(pl_S,ylabel=L"S",yformatter=:auto)
        plot!(pl_L,ylabel=L"Λ",yformatter=:auto)
        plot!(pl_var,ylabel=L"\sigma^2",yformatter=:auto)
        plot!(pl_ac,ylabel=L"ACF(1)",yformatter=:auto)


        #annotations
        annotate!(pl_rr,subfigure_annotation_pos_two_col, text("(a)", :left, annotation_fontsize))
        annotate!(pl_S,subfigure_annotation_pos_two_col, text("(b)", :left, annotation_fontsize))
        annotate!(pl_L,subfigure_annotation_pos_two_col, text("(c)", :left, annotation_fontsize))
        annotate!(pl_var,subfigure_annotation_pos_two_col, text("(d)", :left, annotation_fontsize))
        annotate!(pl_ac,subfigure_annotation_pos_two_col, text("(e)", :left, annotation_fontsize))
    end

    #copy prepared plots for OP
    pl_S_OP = deepcopy(pl_S)
    pl_L_OP = deepcopy(pl_L)
    pl_var_OP = deepcopy(pl_var)
    pl_ac_OP = deepcopy(pl_ac) 

    #------------------------------------read rr intervals---------------------------------
    rrs_data = readdlm(rrs_dir*"rr_$sample.txt")
    rrs_intervals = Float64.(rrs_data[2:end,2]) 
    plot!(pl_rr,rrs_intervals,st=:scatter,lc=linecolor,mc=linecolor,markerstrokewidth=0.0,markershape=:circle,ms=1,ma=0.6);

    #--------------------------------------read and plot measures-----------------------------------------

    #grid
    measures_file_grid = result_files[findall(f -> occursin("measures_grid_$sample", f),result_files)][1]
    @show measures_file_grid
    M_grid = readdlm(data_dir*measures_file_grid)
    window_ends = M_grid[:,1]
    #OP
    measures_file_OP = result_files[findall(f -> occursin("measures_OP_$sample", f),result_files)][1]
    @show measures_file_OP
    M_OP = readdlm(data_dir*measures_file_OP)

    #--------------------------read and search onset times----------------------------
    @show samples[i]
    title!(pl_rr,"$(samples[i])")
    sig_change_data = readdlm(rrs_dir*"ann_sigch_$(samples[i]).txt",skipstart=1)
    rrs_data = readdlm(rrs_dir*"rr_$(samples[i]).txt")
    times = rrs_data[:,1]

    rr_idxs_onset,ann_times_onset = find_annotations(ann_onset,sig_change_data,times)
    rr_idxs_termination,ann_times_termination = find_annotations(ann_term,sig_change_data,times)

    #set plotting interval
    #rr time series indexes start from 1
    #measures indexes start from window_size
    #rr_idxs_onset is the index in the time series vector
    interval_end = min(length(window_ends), rr_idxs_onset[1]-window_size+post_onset_offset)
    interval_start = rr_idxs_onset[1]-pre_onset_offset - window_size
    measures_interval = interval_start:interval_end

    @show rr_idxs_onset[1]
    @show measures_interval
    @show window_ends[measures_interval[1]]
    @show window_ends[measures_interval[end]]
    @show window_ends[measures_interval][1]
    @show window_ends[measures_interval][end]

    #plot grid
    plot!(pl_S,window_ends[measures_interval],M_grid[measures_interval,2],lw=curve_lw,lc=linecolor);
    plot!(pl_L,window_ends[measures_interval],M_grid[measures_interval,3],lw=curve_lw,lc=linecolor);
    plot!(pl_var,window_ends[measures_interval],M_grid[measures_interval,4],lw=curve_lw,lc=linecolor);
    plot!(pl_ac,window_ends[measures_interval],M_grid[measures_interval,5],lw=curve_lw,lc=linecolor); 


    pls = [pl_rr,pl_S,pl_L,pl_var,pl_ac]
    #set xaxis limits
    xlims = (window_ends[measures_interval[1]]-window_size,window_ends[measures_interval[end]])

    #onset line (red) and termination line (green) on subplots
    for pl in pls
        vline!(pl,rr_idxs_onset,ls=:dash,lc=:red,lw=1,label="")
        vline!(pl,rr_idxs_termination,ls=:dash,lc=:green,lw=1,label="")

        #adjust x limits
        #additional offset with window_size on time series plot
        plot!(pl,xlims=xlims,xticks=xticks[i])
    end

    pl = plot(pls...,layout=(5,1),dpi=300)
    push!(plots_grid,pl)

    #plot OP
    plot!(pl_S_OP,window_ends[measures_interval],M_OP[measures_interval,2],lw=curve_lw,lc=linecolor);
    plot!(pl_L_OP,window_ends[measures_interval],M_OP[measures_interval,3],lw=curve_lw,lc=linecolor);
    plot!(pl_var_OP,window_ends[measures_interval],M_OP[measures_interval,4],lw=curve_lw,lc=linecolor);
    plot!(pl_ac_OP,window_ends[measures_interval],M_OP[measures_interval,5],lw=curve_lw,lc=linecolor); 

    pls_OP = [pl_rr,pl_S_OP,pl_L_OP,pl_var_OP,pl_ac_OP]
    #onset line (red) and termination line (green) on subplots
    for pl in pls_OP
        vline!(pl,rr_idxs_onset,ls=:dash,lc=:red,lw=1,label="")
        vline!(pl,rr_idxs_termination,ls=:dash,lc=:green,lw=1,label="")

        #adjust x limits
        #additional offset with window_size on time series plot
        plot!(pl,xlims=xlims,xticks=xticks[i])
    end

    pl = plot(pls_OP...,layout=(5,1),dpi=300)
    push!(plots_OP,pl)

end

fig_dir_name = "figs"
fig_dir = "../../../../" * fig_dir_name * "/" 
!(fig_dir_name in readdir("../../../../")) && (mkdir(fig_dir))

pl = plot(plots_grid...,layout = (1,length(plots_grid)),size=(colfig_size[1]*1.5,colfig_size[2]))
savefig(pl,fig_dir*"fib_selection_paper_sm" * "_grid_$(Int(grid_size))" * "_window_$(Int(window_size))"*".png")

pl_OP = plot(plots_OP...,layout = (1,length(plots_OP)),size=(colfig_size[1]*1.5,colfig_size[2]))
savefig(pl_OP,fig_dir*"fib_selection_paper_sm" * "_OP_$(Int(w))" * "_window_$(Int(window_size))"*".png")



#=
#zoom around onset
for i in 1:length(plots)
    #vertical plots
    for j in 1:length(plots[i])
        xlims!(plots[i][j],(time_spans[i][end]-1000,Inf))
    end
end
pl3_zoom = plot(plots...,layout = (1,length(plots)),size=(2000,1000),leftmargin=12Plots.mm,guidefontsize=18,dpi=300)
plot!(pl3_zoom,plot_title="Window: $window_size")
savefig(pl3_zoom,"figs/LTAFDatabase/summary/fib_selection_paper_zoom_sm"*"_window_$window_size"*".png")
=#