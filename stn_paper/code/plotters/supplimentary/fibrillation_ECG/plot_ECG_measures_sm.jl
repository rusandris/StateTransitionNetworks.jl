using DelimitedFiles
#using CSV,DataFrames
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
grid_sizes = Int.(method_params[1,:])
ws = Int.(method_params[2,:])
window_size,grid_min,grid_max,grid_min06,grid_max06 = method_params[3:7]
grid_edges = [grid_min,grid_max]
grid_edges06 = [grid_min06,grid_max06]
rrs_dir = "../../../../data/supplimentary/fibrillation_ECG/Long_Term_AF_Database/"
#--------------------plot data----------------------

fig_dir = "../../../../figs/" 
mkpath(fig_dir)

subfigure_annotation_pos = (-0.45,1.0)

plot_params = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize-2,
tickfontsize=tickfontsize-5,
titlefontsize=20,
left_margin=reduced_left_margin,
right_margin=reduced_right_margin,
#legend=:none,
xformatter=:none,
yformatter=:auto,
yguidefontrotation=0,
marker_size=marker_size_colored,
dpi=300)

#params for different grid sizes, pattern lengths
linecolor1 = :gray10 #[:gray70,:gray40,:gray10]
linecolor2 = :red
linecolors = [:gray70,:gray40,:gray10]
la = [0.3,0.6,1.0]
#la = [0.3,0.6,1.0][end] #slect only highest line alpha 

#plots are zoomed in on onset
pre_onset_offset = 500 #xaxis limits: nr of idxs before onset
post_onset_offset = 500 #xaxis limits: nr of idxs after onset
post_onset_offset_08 = 200 #xaxis limits: nr of idxs after onset
pre_onset_offset_08 = 200 #xaxis limits: nr of idxs after onset

post_offsets = [post_onset_offset_08,fill(post_onset_offset,4)...]
pre_offsets = [pre_onset_offset_08,fill(pre_onset_offset,4)...]


ann_onset = "(AFIB"
ann_term = "(N"

xticks = [[1600:200:3000;],[4500:200:6500;],[1200:200:2800;],[1400:200:2800;]]
#xticks = [[1200:300:2300;],[1600:300:2800;]]

ylims_rr = (0.0,1.6)
plots_grid = []
plots_OP = []
for (i,sample) in enumerate(samples)

    pre_onset_offset = pre_offsets[i]
    post_onset_offset = post_offsets[i]


    #set rr time series ylims and others
    if i == length(samples)
        ylims_local = grid_edges06
    else
        ylims_local = grid_edges
    end
    #init plots with standard formatting
    pl_rr = plot(;legend=false,ylims=ylims_rr,yticks=[ylims_rr[1]+0.2:0.6:ylims_rr[2]-0.2;],plot_params...,bottom_margin=-2Plots.mm,yguidefontsize=guidefontsize-10)
    pl_S = plot(;legend=false,plot_params...,top_margin=-2Plots.mm,bottom_margin=-2Plots.mm,yguidefontsize=guidefontsize-8) #ylims=(-0.1,2.1),yticks=[0.0,1.0,2.0]
    pl_L = plot(;legend=false,plot_params...,top_margin=-2Plots.mm,bottom_margin=-2Plots.mm) #ylims=(-0.1,3.0),yticks=[0:1:2;]
    pl_var = plot(;legend=false,plot_params...,top_margin=-2Plots.mm,bottom_margin=-2Plots.mm) #ylims=(0.0,0.06),yticks=[0.02,0.03,0.04]
    pl_ac = plot(;legend=false,plot_params...,xformatter=:auto,top_margin=-2Plots.mm,bottom_margin=-2Plots.mm,yguidefontsize=guidefontsize-10) #ylims=(-0.6,1.1),yticks=[-0.4,0.0,0.4]
    xlabel!(pl_ac,"index",xguidefontsize=guidefontsize-12)

    #ylabels and annotation only for the first column 
    if i == 1
        plot!(pl_rr,ylabel=L"\Delta t_{RR}(s)",yformatter=:auto)
        plot!(pl_S,ylabel=L"S,PE",yformatter=:auto)
        plot!(pl_L,ylabel=L"Λ",yformatter=:auto)
        plot!(pl_var,ylabel=L"\sigma^2",yformatter=:auto)
        plot!(pl_ac,ylabel=L"ACF(1)",yformatter=:auto)


        #annotations
        annotate!(pl_rr,subfigure_annotation_pos, text("(a)", :left, annotation_fontsize))
        annotate!(pl_S,subfigure_annotation_pos, text("(b)", :left, annotation_fontsize))
        annotate!(pl_L,subfigure_annotation_pos, text("(c)", :left, annotation_fontsize))
        annotate!(pl_var,subfigure_annotation_pos, text("(d)", :left, annotation_fontsize))
        annotate!(pl_ac,subfigure_annotation_pos, text("(e)", :left, annotation_fontsize))
    end

    #------------------------------------read rr intervals---------------------------------
    rrs_data = readdlm(rrs_dir*"rr_$sample.txt")
    rrs_intervals = Float64.(rrs_data[2:end,2]) 
    plot!(pl_rr,rrs_intervals,st=:scatter,lc=linecolor,mc=linecolor,markerstrokewidth=0.0,markershape=:circle,ms=1.5,ma=0.8);

    #--------------------------------------read and plot measures-----------------------------------------
    #--------------------------read and search onset times----------------------------
    @show samples[i]
    title!(pl_rr,"$(samples[i])")
    sig_change_data = readdlm(rrs_dir*"ann_sigch_$(samples[i]).txt",skipstart=1)
    rrs_data = readdlm(rrs_dir*"rr_$(samples[i]).txt")
    times = rrs_data[:,1] #actual time points in seconds

    #find signal change annotations (fibrillation onset and termination flags)
    rr_idxs_onset,ann_times_onset = find_annotations(ann_onset,sig_change_data,times)
    rr_idxs_termination,ann_times_termination = find_annotations(ann_term,sig_change_data,times)
    @show rr_idxs_onset,ann_times_onset
    @show rr_idxs_termination,ann_times_termination


    #read window ends
    measures_file_grid = result_files[findall(f -> occursin("measures_grid_$sample", f),result_files)][1]
    #use window ends of the first file of the sample
    M_grid = readdlm(data_dir*measures_file_grid)
    window_ends = M_grid[:,1]

    #heartbeat indices starting from 1
    idxs = 1:window_ends[end]

    #set plotting interval
    #rr time series indexes start from 1
    #measures indexes start from window_size
    #rr_idxs_onset is the index in the time series vector

    #loop through discretization params
    for i in 1:length(grid_sizes)
        w = ws[i]
        grid_size = grid_sizes[i] 

        #grid
        measures_file_grid = result_files[findall(f -> occursin("measures_grid_$sample"*"_grid_$grid_size", f),result_files)][1]
        @show measures_file_grid
        M_grid = readdlm(data_dir*measures_file_grid)
        #window_ends = M_grid[:,1]
        #OP
        measures_file_OP = result_files[findall(f -> occursin("measures_OP_$sample"*"_OP_$w", f),result_files)][1]
        @show measures_file_OP
        M_OP = readdlm(data_dir*measures_file_OP)

        @show length(idxs),length(M_grid[:,2])
        #plot grid measures
        plot!(pl_S,window_ends,M_grid[:,2],lw=curve_lw,
            label=L"S(n=%$grid_size)",lc=linecolors[i],la=la[i]); #S
        #plot!(pl_S,window_ends[measures_interval],M_grid[measures_interval,4] ./ log(factorial(w)),lw=curve_lw,
        #    label=L"PE(w=%$w)",lc=linecolor2,ls=:dash,la=la[i]); #C1 
        
        plot!(pl_S,window_ends,M_OP[:,4] ./ log(factorial(w)),lw=curve_lw,
            label=L"PE(w=%$w)",lc=linecolor2,la=la[i]); #C1 OP

        plot!(pl_L,window_ends,M_grid[:,3],lw=curve_lw,
            label=L"\Lambda(n=%$grid_size)",lc=linecolors[i],la=la[i]); #L
        plot!(pl_L,window_ends,M_grid[:,5],lw=curve_lw,
            label=L"C_2(n=%$grid_size)",lc=:purple,la=la[i]); #C2
        plot!(pl_var,window_ends,M_grid[:,6],lw=curve_lw,lc=linecolors[i]); #var
        plot!(pl_ac,window_ends,M_grid[:,7],lw=curve_lw,lc=linecolors[i]); #acf

    end

    #set xaxis limits
    #xlims = (window_ends[measures_interval[1]]-window_size,window_ends[measures_interval[end]])
    #@show xlims
    #onset line (red) and termination line (green) on subplots
    plot!(pl_rr,xlims=(rr_idxs_onset[1]-pre_onset_offset,rr_idxs_onset[1]+post_onset_offset),xticks=xticks[i])
    vline!(pl_rr,rr_idxs_onset,ls=:dash,lc=:red,lw=2,label="")
    vline!(pl_rr,rr_idxs_termination,ls=:dash,lc=:green,lw=2,label="")

    pls = [pl_S,pl_L,pl_var,pl_ac]
    for pl in pls
        
        vline!(pl,rr_idxs_onset,ls=:dash,lc=:red,lw=2,label="")
        vline!(pl,rr_idxs_termination,ls=:dash,lc=:green,lw=2,label="")

        #adjust x limits
        #additional offset with window_size on time series plot
        #plot!(pl,xlims=(1,window_ends[end]),xticks=xticks[i]) #full x scale
        plot!(pl,xlims=(rr_idxs_onset[1]-pre_onset_offset,rr_idxs_onset[1]+post_onset_offset),xticks=xticks[i])
    end
    pl = plot(pl_rr,pls...,layout=(5,1),dpi=300)
    push!(plots_grid,pl)

end


pl = plot(plots_grid...,layout = (1,length(plots_grid)),size=(colfig_size[1]*1.5,colfig_size[2]-300),left_margin=7Plots.mm);
savefig(pl,fig_dir*"ECG_measures_sm" * "_window_$(Int(window_size))"*".pdf")
