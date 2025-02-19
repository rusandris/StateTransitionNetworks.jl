using DelimitedFiles
using CSV,DataFrames
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting_params.jl")

#-----------------------params,constants---------------------------
data_dir = "../../../../data/supplimentary/noisey_timeseries/"
noise_levels = [0.0,0.05,0.1]
ws = [3,5,7] #OP pattern lengths
grid_sizes = [10,32] #grid sizes
xticks_henon = [1.2:0.05:1.4;]

#---------------------henon stuff-------------------
henon_dir = data_dir*"henon_data/"

#grid
measures_dir_grid10 = henon_dir*"henon_entropy_lambda/noise_over_data/henon_noise_over_data_grid_10/"
measures_dir_grid32 = henon_dir*"henon_entropy_lambda/noise_over_data/henon_noise_over_data_grid_32/"

measure_dirs_grid = [measures_dir_grid10,measures_dir_grid32]

filename_template_entropies = "henon_entropies_T1E+06_Ttr1E+03_b_0.3"
filename_template_lambdas = "henon_lambdas_T1E+06_Ttr1E+03_b_0.3"


filenames_entropies_grid = []
filenames_lambdas_grid = []
for grid_size in grid_sizes

    filename_endings = ["_grid_$(grid_size)param_1.2_1.4_noise_0",
        "_grid_$(grid_size)param_1.2_1.4_noise_0.05",
        "_grid_$(grid_size)param_1.2_1.4_noise_0.1" ].* ".txt"

    henon_grid_entropies = filename_template_entropies .*  filename_endings
    push!(filenames_entropies_grid,henon_grid_entropies)

    henon_grid_lambdas = filename_template_lambdas .*  filename_endings
    push!(filenames_lambdas_grid,henon_grid_lambdas)

end

#OP
measures_dir_OP = henon_dir*"henon_ordinal/noise_over_data/henon_noise_over_data_ordinal_"
measures_dirs_OP = measures_dir_OP .* ["3" * "/","5" * "/","7" * "/"]

filename_template_entropies = "henon_entropies_T1E+06_Ttr1E+03_b_0.3_grid_32param_1.2_1.4"
filename_template_lambdas = "henon_lambdas_T1E+06_Ttr1E+03_b_0.3_grid_32param_1.2_1.4"

filenames_entropies_OP = []
filenames_lambdas_OP = []

for w in ws
    filename_endings = ["_noise_0","_noise_0.05","_noise_0.1" ] .* "_ordinal_$w.txt"
    henon_OP_entropies = filename_template_entropies .*  filename_endings
    push!(filenames_entropies_OP,henon_OP_entropies)

    henon_OP_lambdas = filename_template_lambdas .*  filename_endings
    push!(filenames_lambdas_OP,henon_OP_lambdas)
end



#-------------------plots--------------------------
plot_params = (
guidefontsize=guidefontsize,
legendfontsize=legendfontsize,
tickfontsize=tickfontsize,
titlefontsize=20,
left_margin=left_margin,
top_margin=reduced_top_margin,
right_margin=reduced_right_margin,
legend=false,
xformatter=:none,
yformatter=:none,
yguidefontrotation=0,
marker_size=marker_size_colored,
dpi=300)

linealphas = [0.6,1.0]
linealphas_OP = [0.3,0.6,1.0]
rs = readdlm(measure_dirs_grid[1]*filenames_entropies_grid[1][1])[:,1]
S_labels = ["(a)","(b)","(c)"]
S_OP_labels = ["(d)","(e)","(f)"]
L_labels = ["(g)","(h)","(i)"]
L_OP_labels = ["(j)","(k)","(l)"]

plots_henon_S = []
plots_henon_L = []
plots_henon_S_OP = []
plots_henon_L_OP = []
for (i,noise_level) in enumerate(noise_levels)

    local_annot_pos = subfigure_annotation_pos_one_col

    #skip y labels
    if i == 1
        ylabel_S_local=L"S"
        ylabel_S_local_OP=L"S_{\mathrm{OP}}"
        ylabel_L_local=L"\Lambda"
        ylabel_L_local_OP=L"\Lambda_{\mathrm{OP}}" 
    else
        ylabel_S_local=""
        ylabel_L_local=""
        ylabel_S_local_OP=""
        ylabel_L_local_OP=""
    end

    #grid
    pl_S = plot(;ylims=(0.0,2.0),yticks=[0:2;],title = L"\sigma_{\mathrm{n}}/\sigma_{\mathrm{s}} = %$noise_level",xticks=xticks_henon, ylabel=ylabel_S_local,plot_params...)
    pl_L = plot(;ylims=(0.0,1.5),yticks=[0:0.5:1.5;],title = L"\sigma_{\mathrm{n}}/\sigma_{\mathrm{s}} = %$noise_level",xticks=xticks_henon,ylabel=ylabel_L_local,plot_params...)

    #OP
    pl_S_OP = deepcopy(pl_S)
    pl_L_OP = deepcopy(pl_L)

    #re-enable yticks
    if i == 1
        plot!(pl_S,yformatter=:auto,legend=true)
        plot!(pl_L,yformatter=:auto,legend=true)
        plot!(pl_S_OP,yformatter=:auto,legend=true)
        plot!(pl_L_OP,yformatter=:auto,legend=true)
        local_annot_pos = subfigure_annotation_pos_two_col_alter_layout
    end

    annotate!(pl_S,local_annot_pos, text(S_labels[i], :left, annotation_fontsize))
    annotate!(pl_L,local_annot_pos, text(L_labels[i], :left, annotation_fontsize))
    annotate!(pl_S_OP,local_annot_pos, text(S_OP_labels[i], :left, annotation_fontsize))
    annotate!(pl_L_OP,local_annot_pos, text(L_OP_labels[i], :left, annotation_fontsize))

    plot!(pl_S_OP;ylabel=ylabel_S_local_OP,ylims=(0.0,1.0),yticks=[0:0.5:1.0;])
    plot!(pl_L_OP;ylabel=ylabel_L_local_OP,ylims=(0.0,2.0),yticks=[0:2;],xformatter=:auto,xlabel=L"a")

    #grid
    for (j,grid_size) in enumerate(grid_sizes)
        measures_dir = measure_dirs_grid[j]
        @show measures_dir
        Ss = readdlm(measures_dir*filenames_entropies_grid[j][i])[:,2]
        Λs = readdlm(measures_dir*filenames_lambdas_grid[j][i])[:,2]
    
        plot!(pl_S,rs,Ss,lc=:gray10,label=L"n = %$grid_size"*L"^2",la=linealphas[j])
        plot!(pl_L,rs,Λs,lc=:gray10,label=L"n = %$grid_size"*L"^2",la=linealphas[j])
    end
    push!(plots_henon_S,pl_S)
    push!(plots_henon_L,pl_L)


    #OP
    for (j,w) in enumerate(ws)
        measures_dir = measures_dirs_OP[j]
        @show measures_dir
        Ss = readdlm(measures_dir*filenames_entropies_OP[j][i])[:,2]
        Λs = readdlm(measures_dir*filenames_lambdas_OP[j][i])[:,2]
    
        plot!(pl_S_OP,rs,Ss,lc=:gray10,label=L"w = %$w",la=linealphas_OP[j])
        plot!(pl_L_OP,rs,Λs,lc=:gray10,label=L"w = %$w",la=linealphas_OP[j])
    end
    push!(plots_henon_S_OP,pl_S_OP)
    push!(plots_henon_L_OP,pl_L_OP)

end

pl_henon_S = plot(plots_henon_S...,plots_henon_S_OP...,layout=(2,3))
pl_henon_L = plot(plots_henon_L...,plots_henon_L_OP...,layout=(2,3))
pl_henon = plot(pl_henon_S,pl_henon_L,layout=(2,1),size=bigfig_size)

