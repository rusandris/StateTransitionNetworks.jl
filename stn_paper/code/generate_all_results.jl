#----------------------reproduce all results--------------------------
#some parameter might be set to quick testing 
cd(@__DIR__)

include("functions_utils.jl")

#main results 
include("logistic_map_main.jl")
include("henon_map_main.jl")
#save inset data separately for main
include("save_inset_data_logistic.jl")
include("save_inset_data_henon.jl")

#supplimentary
include("logistic_map_sm.jl")
include("henon_map_sm.jl")
include("calc_save_var_acf_sm.jl")
include("nonauto_logistic_sm.jl")
include("ECG_data_measures_sm.jl")

