Our package `StateTransitionNetworks.jl` and other functions and scripts to generate figures are written in the Julia programming language.

For Julia installation instructions see: https://julialang.org/downloads/

Run everything from the `stn_paper/` dir as working directory to avoid path errors.

#### INSTALL AND SETUP

1. Activate the environment:
    * in the `stn_paper/` directory open a Julia `REPL`, activate the env with: `] activate .` (the `]`is needed to enter package mode (blue) of the `REPL` and only needed once)

2. Add `StateTranstionNetworks.jl` to the env:

   * if you downloaded the whole repository: ` add ../../StateTransitionNetworks.jl/`

   * if you only have the `stn_paper/` dir with the scripts, or you use Windows: 

     ` add https://github.com/rusandris/StateTransitionNetworks.jl.git`

3. Instantiate to install all packages: ` instantiate`

#### RUN

1. Run from `REPL`: `include(code/code.jl)`
2. Run in background (Linux): `julia --project=./ code/code.jl < /dev/null > logfile.log 2>&1 &`
3. Run from `VSCode`: select environment and run

Note: The scripts save the values continuously in `data/` and log status at every 100th parameter value in the stdout/logfiles. If you run the same script multiple times, they all save data in the same files (duplicates can occur, and this can cause errors in plotting)

#### PLOT

To get the main figure, run `code/plotters/plot_maps_main.jl`. See the output in `figs/`

It is advised to do this in a separate Julia session, if you intend to run the plotter from the `REPL`.
