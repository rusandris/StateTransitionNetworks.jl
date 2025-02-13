#---------------------global plotting params---------------------
global guidefontsize::Int = 25
global legendfontsize::Int = 14
global tickfontsize::Int = 14
global annotation_fontsize::Int = 24 
global reduced_right_margin = 2Plots.mm
global reduced_left_margin = 4Plots.mm
global left_margin = 9Plots.mm
global reduced_top_margin = -6Plots.mm
global reduced_top_margin_sm = -2Plots.mm
global top_margin = 5Plots.mm
global halfcolfig_size = (500,800) #not sure of it matters
global colfig_size = (1000,1200)
global bigfig_size = (1500,1200)
global marker_size_colored = 12 
global ylims_od = (-0.2,1.45)
global inset_box = bbox(0.63,-0.05,0.2,0.3)
global subfigure_annotation_pos = (-0.35,1.0)
global subfigure_annotation_pos_henon = (-0.13,1.0)
global subfigure_annotation_pos_sm = (-0.15,1.0)
global subfigure_annotation_pos_sm_henon = (-0.1,1.0)
global subfigure_annotation_pos_alter = (-0.3,1.0)
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
global h_top_henon::Float64 = 0.465
global λ_henon::Float64 = 0.42
global ms_od_henon = 0.5
global ma_od_henon = 0.9
global ms_od_logistic = 0.3
global ma_od_logistic = 0.8

global ensemble::Int64 = 10^3 
global N_steps::Int64 = 500
gauss(x,μ,σ) = 1/(σ*sqrt(2π)) * exp(-(x - μ)^2/(2σ^2)) 