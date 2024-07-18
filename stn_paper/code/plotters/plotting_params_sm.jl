#---------------------global plotting params---------------------
global guidefontsize::Int = 25
global legendfontsize::Int = 18
global tickfontsize::Int = 18
global annotation_fontsize::Int = 24 
global reduced_right_margin = 2Plots.mm
global reduced_left_margin = 4Plots.mm
global left_margin = 9Plots.mm
global marker_size_colored = 12 
global subfigure_annotation_pos = (-0.2,0.95)
global subfigure_annotation_pos_henon = (-0.08,0.95)
global inset_tickfontsize = 14
global legend_position = (0.6,0.7)
global marker_shapes = [:rect,:star5,:circle,:diamond]
global marker_colors = [:red,:blue,:green,:orange]
global wl_colors = [:red,:green,:orange]
global alphas = [0.4:0.2:1.0;]
global marker_pos_henon = [0.005,0.02,0.02]
global marker_pos_log = [0.007,0.02,0.02]
global histogram_alpha = 0.5
global nr_bins = 40
global common_y_lims = (0.0,0.075)
global curve_lw = 2

global ensemble::Int64 = 10^3 
global N_steps::Int64 = 500
gauss(x,μ,σ) = 1/(σ*sqrt(2π)) * exp(-(x - μ)^2/(2σ^2)) 