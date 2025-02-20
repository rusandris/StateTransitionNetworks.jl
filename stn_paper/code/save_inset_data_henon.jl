using DelimitedFiles
cd(@__DIR__)
data_dir = "../data_test/main/henon_data/"
inset_param_range = [1.2,1.25] #inset params
result_files = readdir(data_dir)

entropies_file = result_files[findall(f -> occursin("entrop", f),result_files)][1]
lambdas_file = result_files[findall(f -> occursin("lambda", f),result_files)][1]

#cut and save measures to put on the inset
Ss_data = readdlm(data_dir * entropies_file)[:,1:end]
Λs_data = readdlm(data_dir * lambdas_file)[:,1:end]
ps = Ss_data[:,1] #params

#find idxs of inset range
start_idx = findall(p -> p == inset_param_range[1],ps)[1]
stop_idx = findall(p -> p == inset_param_range[2],ps)[1]

ps_inset = ps[start_idx:stop_idx]
Ss = Ss_data[start_idx:stop_idx,5] #analytic entropy (highest order)
Λs = Λs_data[start_idx:stop_idx,5] #analytic lambda (highest order)
Ss_rw = Ss[start_idx:stop_idx,end] #random walk entropy (highest order)
Λs_rw = Λs_data[start_idx:stop_idx,end] #random walk lambda (highest order)

writedlm(data_dir*"henon_inset_measures.txt",hcat(ps_inset,Ss,Ss_rw,Λs,Λs_rw))

