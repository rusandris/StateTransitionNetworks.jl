Our package `StateTransitionNetworks.jl` and other functions and scripts to generate figures are written in the Julia programming language.

These examples are likely but not guaranteed to work with the latest version of the package. **To make sure they work, use the `v0.4.2-pub` release.**

For Julia installation instructions see: https://julialang.org/downloads/

#### INSTALL AND SETUP

Once you have successfully downloaded the `stn_paper/` directory (or by clicking on the `Download ZIP` option on Github, or by cloning the repository `git clone https://github.com/rusandris/StateTransitionNetworks.jl.git --branch v0.4.2-pub`  ), navigate to `stn_paper/` directory and follow these steps:

1. Activate the environment:
    * in the `stn_paper/` directory open a Julia `REPL`, activate the env with: `] activate .` (the `]`is needed to enter package mode (blue) of the `REPL` and only needed once)

2. Add `StateTranstionNetworks.jl` to the env (if you haven't done before into a global Julia env):

   * if you only have the `stn_paper/` dir with the scripts, or you use Windows: 

     ` add https://github.com/rusandris/StateTransitionNetworks.jl.git #v0.4.2-pub  `

   * if you downloaded the whole repository: ` add ../../StateTransitionNetworks.jl/`

     Note: Adding our package is done like this because`StateTransitionNetworks.jl` isn't in the official Julia registry yet.

3. Instantiate to install all packages: ` instantiate`

**All scripts that generate numerical results read and save results to `stn_paper/data`, and plotters save figures to `stn_paper/figs`. **  

#### Generate figures from pre-saved results

Pre-saved results are in `presaved_results.zip` . 

1. Extract to `data/` directory to `stn_paper` (extract as folder `data`)
2. run `code/plotters/plot_all_figs.jl`

The figures from the paper (all except Fig.1 in the Supplimental Material) can be reproduced by running `code/plotters/plot_all_figs.jl` and they will be saved to `stn_paper/figs`.

Figures `logistic_noise_results.png` and `henon_noise_results.png` can only be reproduced using the presaved data at the moment, the source code that generate it can be made available later.

Source code to critical and tent map results (Fig.1) can be made available on request.  

#### Generate results

1. extract `data/` dir (extract as folder `data`) from `data_LTAF_misc.zip` that contains data taken from the LTAF dataset (https://physionet.org/content/ltafdb/1.0.0/). 

2. open a new Julia session

3. To generate the results run source code:

* `generate_all_results.jl`

and plot with

* `plot_all_figs.jl`.

Note that each script changes the working directory to their own, so keep that in mind while running things in the REPL. For example when starting with `stn_paper` as the working directory:

```include("code/generate_all_results.jl")```
```include("plotters/plot_all_figs.jl")```

*Note:* Keep in mind that some method parameters (for example length of calculated time series) are chosen to produce results quickly and of course may lead to less accurate calculation. 

To reproduce the exact figures,  one should use the parameters values in the comments (there are some values that are commented out after certain parameters) in scripts in `/code`:

* `logistic_map_main.jl`
* `henon_map_main.jl`
* `logistic_map_sm.jl`
* `henon_map_sm.jl`
* `calc_save_var_acf_sm.jl`

To run any file from the `code/` directory, here are some options:

* Run from `REPL`: `include("code/logistic_map_main.jl")`

* Run in background (Linux): `julia --project=./ code/logistic_map_main.jl < /dev/null > logfile.log 2>&1 &`

* Run from `VSCode`: select environment and run

Note:  The scripts 

* `logistic_map_main.jl`
* `logistic_map_sm.jl`
* `henon_map_main.jl`
* `henon_map_sm.jl`

save the values continuously in `data/` and log status in the stdout/logfiles. If you run the same source code multiple times in a row, results get appended to the same file, leading to errors.

