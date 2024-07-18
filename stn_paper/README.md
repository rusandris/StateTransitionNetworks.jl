Our package `StateTransitionNetworks.jl` and other functions and scripts to generate figures are written in the Julia programming language.

For Julia installation instructions see: https://julialang.org/downloads/

Run everything from the `stn_paper/` dir as working directory to avoid path errors.

#### INSTALL AND SETUP

Once you have successfully downloaded the `stn_paper/` directory (or by clicking on the `Download ZIP` option on Github, or by cloning the repository `git clone https://github.com/rusandris/StateTransitionNetworks.jl.git `), navigate to `stn_paper/` directory and follow these steps:

1. Activate the environment:
    * in the `stn_paper/` directory open a Julia `REPL`, activate the env with: `] activate .` (the `]`is needed to enter package mode (blue) of the `REPL` and only needed once)

2. Add `StateTranstionNetworks.jl` to the env (if you haven't done before into a global Julia env):

   * if you only have the `stn_paper/` dir with the scripts, or you use Windows: 

     ` add https://github.com/rusandris/StateTransitionNetworks.jl.git`

   * if you downloaded the whole repository: ` add ../../StateTransitionNetworks.jl/`

     Note: Adding our package is done like this because`StateTransitionNetworks.jl` isn't in the official Julia registry yet.

3. Instantiate to install all packages: ` instantiate`

#### RUN

To run any file from the `code/` directory, here are some options:

* Run from `REPL`: `include(code/logistic_map_main.jl)`

* Run in background (Linux): `julia --project=./ code/logistic_map_main.jl < /dev/null > logfile.log 2>&1 &`

* Run from `VSCode`: select environment and run

To generate the results run the following source codes:

* `logistic_map_main.jl`
* `henon_map_main.jl`

Note: To reproduce the exact figures,  one should use the parameters values in the comments (there are some values that are commented out after certain parameters). The ones that are currently set are chosen to produce results quickly and of course may lead to less accurate calculation. . 

Note:  The scripts save the values continuously in `data/` and log status in the stdout/logfiles. If you run the same script multiple times, they all save data in the same files (duplicates can occur, and this can cause errors in plotting).

#### PLOT

To get the main figure, run `code/plotters/plot_maps_main.jl`. 

To get the supplementary figure, run `code/plotters/plot_maps_sm.jl`.

See the output in `figs/` (png).

It is advised to do this in a separate Julia session, if you intend to run the plotter from the `REPL`.
