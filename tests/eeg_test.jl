using StateTransitionNetworks

grid = 20
plane = (1,0)
ens = 1000
steps = 10000
filename = "BB0008_SP2_stg2_F3C3P3O1LOCA2EMG_7856640_7987712.dat"

eeg_data = read_bin(filename,Float32,6)
stn,retcode = stn_analysis(eeg_data;grid=grid,plane=plane,idxs=[2,3],return_stn=true)

@show retcode

lyap, entr = network_measures(stn,ens,steps)
@show lyap, entr

fn, ext = splitext(filename)
plot_stn(stn;filename=fn*".pdf",linetype="curve")

