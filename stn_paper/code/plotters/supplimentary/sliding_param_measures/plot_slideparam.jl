using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)
include("../../plotting_params.jl")

plot_params = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
titlefontsize=20,
left_margin=left_margin,
top_margin=reduced_top_margin,
right_margin=reduced_right_margin,
legend=:none,
xformatter=:none,
yguidefontrotation=0,
marker_size=marker_size_colored,
dpi=300)

#-------------------------------------read results----------------------------------
data_dir = "../../../../data/supplimentary/sliding_param_measures/"

result_files = readdir(data_dir)
entropies_file = result_files[findall(f -> occursin("entropies", f),result_files)][1]
lambdas_file = result_files[findall(f -> occursin("lambdas", f),result_files)][1]
vars_file = result_files[findall(f -> occursin("vars", f),result_files)][1]
acs_file = result_files[findall(f -> occursin("acs", f),result_files)][1]
ts_file = result_files[findall(f -> occursin("timeseries", f),result_files)][1]
time_param_file = result_files[findall(f -> occursin("time_param", f),result_files)][1]
method_param_file = result_files[findall(f -> occursin("method_param", f),result_files)][1]
system_param_file = result_files[findall(f -> occursin("system_param", f),result_files)][1]

Ss = readdlm(data_dir*entropies_file)[:,2:end]
Λs = readdlm(data_dir*lambdas_file)[:,2:end]
vars = readdlm(data_dir*vars_file)[:,2:end]
acs = readdlm(data_dir*acs_file)[:,2:end]
t_windows = Int.(readdlm(data_dir*entropies_file)[:,1])
ts = readdlm(data_dir*ts_file)
time_param = readdlm(data_dir*time_param_file)
method_param = readdlm(data_dir*method_param_file)
system_param = readdlm(data_dir*system_param_file)


#-------------------------------------plots--------------------------------------------------

times = ts[:,1]
rs = ts[:,2]
timeseries = ts[:,3]
n_min,n_max,n_crit = round.(Int,time_param)
n_trans = round(Int,n_min + n_crit)
grid_size,window_size = Int.(method_param)
r_min,r_max,r_crit = system_param
nr_ics = size(Ss)[2]
top_xaxis_ticks = 10 .^ [0,1,2,3,4,5]

#measures_xlims = (rs[1],rs[end])
measures_xlims = (times[1],times[end])
#linealphas = range(0.3,1.0;length=nr_ics)
linealphas = fill(0.3,nr_ics-1)
push!(linealphas,1.0)
lws=[1,1,2]
r_windows = rs[end-length(t_windows)+1:end] 

#------------------------------------------------plot timeseries------------------------------------------
pl_ts = plot(;plot_params...,top_margin=4Plots.mm);
twiny(pl_ts) #shared y axis, multiple x axis
plot!(pl_ts[2],xlims=(rs[1],rs[end]+0.001),xlabel=L"r(t)",xticks=[3.8,3.81,3.82,3.83,3.84],widen=true,guidefontsize=guidefontsize,tickfontsize=tickfontsize)

plot!(pl_ts[1],times,timeseries,st=scatter,ms=0.8,ylims=(0.0,1.2),yticks=[0.0:0.5:1.0;],widen=true,ylabel=L"x(t)",mc=:gray10,ma=0.5,markerstrokewidth=0.0,legend=false);
vline!(pl_ts[1],[n_min+n_crit],ls=:dash,lw=2,lc=:red);


pl_S = plot(;plot_params...);
pl_L = plot(;plot_params...,ylims=(0.0,2.5),yticks=[0.0:1.0:2.5;]);
pl_var = plot(;plot_params...,ylims=(0.0,0.12),yticks=[0.0:0.05:0.1;]);
pl_ac = plot(;ylims=(-0.7,-0.45),yticks=[-0.7:0.1:-0.5;],xlims=[times[1],times[end]],xticks=[times[1]:30000:times[end];],xlabel=L"t",plot_params...,xformatter=:auto);


for i in 1:nr_ics
    #lambda
    plot!(pl_L,t_windows,Λs[:,i],lc=:gray10,lw=lws[i],ylabel=L"\Lambda",xlims = measures_xlims,la = linealphas[i],widen=true,legend=false);
    vline!(pl_L,[n_trans],ls=:dash,lw=2,lc=:red); #red line
    plot!(pl_L,[n_trans],[Λs[n_trans - n_min - window_size,nr_ics]],st=:scatter,mc=:gray10,ms=5,markerstrokewidth=0.0)
    hline!(pl_L,[Λs[n_trans - n_min - window_size,nr_ics]],ls=:dash,lc=:gray50)
    #var
    plot!(pl_var,t_windows,vars[:,i],lc=:gray10,lw=lws[i],ylabel=L"\sigma^2",widen=true,xlims = measures_xlims,la = linealphas[i],legend=false);
    vline!(pl_var,[n_trans],ls=:dash,lw=2,lc=:red); #red line
    plot!(pl_var,[n_trans],[vars[n_trans - n_min - window_size,nr_ics]],st=:scatter,mc=:gray10,ms=5,markerstrokewidth=0.0)
    hline!(pl_var,[vars[n_trans - n_min - window_size,nr_ics]],ls=:dash,lc=:gray50)
    #ac
    plot!(pl_ac,t_windows,acs[:,i],lc=:gray10,lw=lws[i],ylabel=L"ACF(1)",widen=true,xlims = measures_xlims,la = linealphas[i],legend=false);
    vline!(pl_ac,[n_trans],ls=:dash,lw=2,lc=:red); #red line
    plot!(pl_ac,[n_trans],[acs[n_trans - n_min - window_size,nr_ics]],st=:scatter,mc=:gray10,ms=5,markerstrokewidth=0.0)
    hline!(pl_ac,[acs[n_trans - n_min - window_size,nr_ics]],ls=:dash,lc=:gray50)
end

annotate!(pl_ts,subfigure_annotation_pos_one_col, text("(a)", :left, annotation_fontsize))
annotate!(pl_L,subfigure_annotation_pos_one_col, text("(b)", :left, annotation_fontsize))
annotate!(pl_var,subfigure_annotation_pos_one_col, text("(c)", :left, annotation_fontsize))
annotate!(pl_ac,subfigure_annotation_pos_one_col, text("(d)", :left, annotation_fontsize))


pl_log = plot(pl_ts,pl_L,pl_var,pl_ac,layout=(4,1),size=colfig_size);
#plot!(plot_title = "Logistic map",legend=false,guidefontsize=16,tickfontsize=12,legendfontsize=12,leftmargin=18Plots.mm,rightmargin=10Plots.mm,dpi=300);

#------------------------------------------save fig----------------------------------------------
fig_dir_name = "figs"
fig_dir = "../../../../" * fig_dir_name * "/" 
!(fig_dir_name in readdir("../../../../")) && (mkdir(fig_dir))
savefig(pl_log,fig_dir*"logistic_map_measures_comparison_slideparam_many_ics_zoom_tslength$n_max" * "_window$window_size" * "_grid$grid_size" * "_r_min$r_min" * "_r_max$r_max" * ".png")
