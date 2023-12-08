# CLiMA-Hybrid Repository

## Repository to run CLiMA-Hybrid model

Just activate the environment in the root folder [CLiMA-Hybrid] and do:

```julia
] add https://github.com/lazarusA/Land.git#single_site
```

or, if you want to do further changes to `Land` locally then do 

```julia
] dev --local https://github.com/lazarusA/Land.git
```

`cd` to `Land` and checkout the branch `single_site`, i.e.

```shell
> cd Land
Land> git checkout single_site
Land> git pull # get the latest changes if any 
```
start working with your scripts by running them or adding more packages to your env `CLiMA-Hybrid`.

# Download files
The examples here need the following files:

```
$ wget https://github.com/CliMA/Land/raw/v0.1/examples/debug.jld2
$ wget https://github.com/CliMA/Land/raw/v0.1/examples/debug.nc
```

into `ClimaLand_examples`.

### To run CLiMA example.jl:
0) Simply go to ClimaLand_examples/example.jl and start running line by line (or all the code) in vscode.
1) Or, initialize the Julia REPL, then activate and instantiate the Julia env and then do:
2) ``` cd ClimaLand_examples```
3)  ``` include("example1.jl") ``` for standard REPL CLiMA run.
4) Merge the output from the standard run to the input file debug.nc, to have a complete file, used for training the hybrid model.

### To run the hybrid model predicting leaf stomatal conductance and transpiration:
1) ``` cd Experiments  ```
2) ``` include("hyb_model_ex3a.jl") ``` for standard hybrid model run.

### To run the hybrid model predicting canopy-level conductance and transpiration:
1) ``` cd Experiments  ```
2) ``` include("hyb_model_ex3a_gc.jl") ``` for standard hybrid model run.

### To run hybrid model in CLiMA to predict at leaf-level:
1) cd  ``` ClimaLand_examples```
2) Run  ``` include("example_hyb2.jl") ```

### To run hybrid model in CLiMA to predict at canopy-level:
1) cd  ``` ClimaLand_examples```
2) Run  ``` include("example_hyb2_gc.jl") ```
